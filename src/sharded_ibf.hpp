#pragma once

// Included after the v2 format-3 archive functions. Kept separate so physical
// storage cannot accidentally introduce shard-local biological decisions.
#include "trusted_archive.hpp"

double elapsed(Clock::time_point start);
std::size_t peak_rss_bytes();
double observed_ibf_union_estimate(Ibf const &, ReferenceInfo const &, Options const &);

void plan_shards(Options & o, ReferenceInfo const & info) {
    auto bits = planned_bin_bits(o.fragment_size,o.hash_functions,o.fpr);
    auto single = planned_index_bytes(info.bin_count,bits);
    if (o.ibf_layout == "auto") {
        o.ibf_layout = info.bin_count > o.max_bins_per_ibf ||
            (o.ibf_max_index_memory_gb > 0 && single > o.ibf_max_index_memory_gb*1e9) ? "sharded" : "single";
        fraglncrna::production_log.event("shard-plan","stage_progress","auto decision=" + o.ibf_layout +
            " single_bytes=" + std::to_string(single) + " real_bin_cap=" + std::to_string(o.max_bins_per_ibf));
    }
    if (o.ibf_build_only && o.ibf_layout == "single" && o.save_index.empty())
        throw std::runtime_error("Auto selected single for build-only; use --ibf-layout sharded or explicitly select single with --save-index");
    auto cap = o.ibf_layout == "sharded" ? o.max_bins_per_ibf : info.bin_count;
    auto count = info.bin_count/cap + (info.bin_count%cap != 0);
    auto full = planned_index_bytes(std::min(cap,info.bin_count),bits);
    auto tail = planned_index_bytes(info.bin_count-(count-1)*cap,bits);
    if (count-1 > (std::numeric_limits<std::size_t>::max()-tail)/full)
        throw std::runtime_error("Shard disk estimate overflow");
    auto total = full*(count-1)+tail;
    auto resident = o.ibf_shard_storage == "memory" ? total : full;
    if (o.ibf_max_index_memory_gb > 0 && resident > o.ibf_max_index_memory_gb*1e9)
        throw std::runtime_error("Planned resident IBF bytes " + std::to_string(resident) +
            " exceed --ibf-max-index-memory-gb; lower --max-bins-per-ibf or use disk storage");
    // Memory mode must not rely on an unlimited default to exhaust the host.
    auto available = double(fraglncrna::memory_limit_bytes());
    if (resident > available*0.5)
        throw std::runtime_error("Planned resident IBF exceeds half physical RAM; lower --max-bins-per-ibf and use sharded disk storage");
    auto destination = o.ibf_layout == "sharded" ? fs::absolute(o.ibf_shard_dir) : fs::absolute(o.output).parent_path();
    while (!fs::exists(destination)) destination = destination.parent_path();
    auto free = fs::space(destination).available;
    if (o.ibf_layout == "sharded" && free < total + full + 1048576)
        throw std::runtime_error("Insufficient free disk for durable shards plus atomic staging; choose --ibf-shard-dir on a larger filesystem");
    fraglncrna::production_log.event("shard-plan","stage_completed","layout=" + o.ibf_layout +
        " contigs=" + std::to_string(info.ids.size()) + " bases=" + std::to_string(info.total_bases) +
        " k=" + std::to_string(o.k) + " fragment_size=" + std::to_string(o.fragment_size) +
        " valid_windows=" + std::to_string(info.valid_windows) + " global_bins=" + std::to_string(info.bin_count) +
        " shards=" + std::to_string(count) + " max_real_bins=" + std::to_string(std::min(cap,info.bin_count)) +
        " max_padded_bins=" + std::to_string(((std::min(cap,info.bin_count)+63)/64)*64) +
        " estimated_disk_bytes=" + std::to_string(total) + " estimated_resident_index_bytes=" + std::to_string(resident) +
        " free_disk_bytes=" + std::to_string(free) + " retention=" + o.ibf_shard_retention +
        " memory_bound_scope=IBF_bitvectors_not_total_process_RSS");
}

void sync_shard_directory(fs::path const & directory) {
    int descriptor = open(directory.c_str(),O_RDONLY | O_DIRECTORY | O_CLOEXEC);
    if (descriptor < 0) throw std::runtime_error("Cannot open shard directory for durable commit: " + directory.string());
    auto result = fsync(descriptor); close(descriptor);
    if (result) throw std::runtime_error("Cannot sync shard directory; select a filesystem supporting durable index commits");
}

// Atomic self-checksummed text records avoid a two-file manifest/checksum commit.
void atomic_shard_text(fs::path const & path, std::string const & text) {
    auto temporary = unique_staging_path(path);
    try {
        std::ofstream out(temporary,std::ios::binary); out << text; out.close();
        if (!out) throw std::runtime_error("Cannot write shard transaction " + temporary.string());
        int fd = open(temporary.c_str(),O_RDONLY | O_NOFOLLOW);
        if (fd < 0) throw std::runtime_error("Cannot open shard transaction for sync");
        int result = fsync(fd); close(fd);
        if (result) throw std::runtime_error("Cannot flush shard transaction");
        fs::rename(temporary,path); sync_shard_directory(path.parent_path());
    } catch (...) { std::error_code ignored; fs::remove(temporary,ignored); throw; }
}

class ShardDirectoryLock {
    int fd_ = -1;
public:
    explicit ShardDirectoryLock(fs::path const & dir) {
        fs::create_directories(dir);
        fd_ = open((dir/".lock").c_str(),O_CREAT | O_RDWR | O_NOFOLLOW,0600);
        if (fd_ < 0) throw std::runtime_error("Cannot open shard lock in " + dir.string());
        if (flock(fd_,LOCK_EX | LOCK_NB)) {
            close(fd_); fd_ = -1;
            throw std::runtime_error("Shard directory locked by another run; wait or choose another --ibf-shard-dir");
        }
    }
    ~ShardDirectoryLock() { if (fd_ >= 0) { flock(fd_,LOCK_UN); close(fd_); } }
    ShardDirectoryLock(ShardDirectoryLock const &) = delete;
    ShardDirectoryLock & operator=(ShardDirectoryLock const &) = delete;
};

class ShardedIbf {
    struct Shard {
        std::size_t begin = 0, bins = 0, bytes = 0;
        std::string file, sha, provenance;
        bool owned = false;
        std::unique_ptr<Ibf> memory;
        std::unique_ptr<TrustedArchive> trusted, local_trusted;
        std::uint64_t last_use = 0;
    };
    Options const & o_;
    ReferenceInfo & info_;
    std::string bed_sha_, prefix_;
    fs::path dir_;
    ShardDirectoryLock lock_;
    fs::path local_root_, local_dir_;
    std::unique_ptr<ShardDirectoryLock> local_lock_;
    std::vector<Shard> shards_;
    std::vector<std::size_t> base_offsets_;
    double guard_bound_ = 0;
    std::size_t cache_budget_ = 0, cache_bytes_ = 0, resident_peak_ = 1;
    std::uint64_t cache_clock_ = 0;
    void make_room(Shard const & x) {
        if (o_.ibf_shard_storage != "disk") return;
        auto bytes = planned_index_bytes(x.bins,info_.bin_bits);
        auto budget = std::max(cache_budget_,bytes);
        while (cache_bytes_ > budget-bytes) {
            Shard * oldest = nullptr;
            for (auto & other : shards_) if (other.memory && (!oldest || other.last_use < oldest->last_use)) oldest = &other;
            if (!oldest) throw std::runtime_error("Internal cache eviction failure");
            cache_bytes_ -= planned_index_bytes(oldest->bins,info_.bin_bits);
            oldest->memory.reset(); ++cache_evictions;
        }
    }
    void retain(Shard & x, std::unique_ptr<Ibf> ibf) {
        auto bytes = planned_index_bytes(x.bins,info_.bin_bits);
        if (o_.ibf_shard_storage == "disk" && bytes > cache_budget_) return;
        if (o_.ibf_shard_storage == "disk") {
            while (cache_bytes_ > cache_budget_-bytes) {
                Shard * oldest = nullptr;
                for (auto & other : shards_) if (other.memory && (!oldest || other.last_use < oldest->last_use)) oldest = &other;
                if (!oldest) throw std::runtime_error("Internal shard cache accounting failure");
                cache_bytes_ -= planned_index_bytes(oldest->bins,info_.bin_bits);
                oldest->memory.reset(); ++cache_evictions;
            }
        }
        x.memory = std::move(ibf); x.last_use = ++cache_clock_; cache_bytes_ += bytes;
        std::size_t count = 0; for (auto const & other : shards_) count += bool(other.memory);
        resident_peak_ = std::max(resident_peak_,count);
        cache_peak_bytes = std::max(cache_peak_bytes,cache_bytes_);
    }
    fs::path read_path(Shard const & x) const {
        return x.local_trusted ? x.local_trusted->pinned_path() : x.trusted->pinned_path();
    }
    void check(Shard const & x) const {
        x.trusted->check(); if (x.local_trusted) x.local_trusted->check();
    }
    std::size_t local_used_bytes() const {
        std::size_t bytes = 0;
        for (auto const & item : fs::recursive_directory_iterator(local_root_)) {
            if (item.is_symlink()) throw std::runtime_error("Symlink in local read-cache root refused");
            if (item.is_regular_file()) {
                auto count = item.file_size();
                if (count > std::numeric_limits<std::size_t>::max()-bytes) throw std::runtime_error("Local cache byte counter overflow");
                bytes += count;
            }
        }
        return bytes;
    }
    void stage_local(Shard & x) {
        if (local_dir_.empty()) return;
        auto destination = local_dir_/(x.sha+".ibf");
        if (fs::is_symlink(destination)) throw std::runtime_error("Symlink local read-cache entry refused");
        if (!fs::exists(destination)) {
            auto used = local_used_bytes();
            auto budget = o_.ibf_local_cache_max_mb*1048576;
            if (x.bytes > budget || used > budget-x.bytes) throw std::runtime_error("Local read-cache byte budget exhausted");
            if (fs::space(local_dir_).available < x.bytes+1048576) throw std::runtime_error("Insufficient disk space for explicit local index cache");
            auto temporary = unique_staging_path(destination);
            try {
                fs::copy_file(x.trusted->pinned_path(),temporary,fs::copy_options::overwrite_existing);
                x.trusted->check(); fs::rename(temporary,destination); local_copy_bytes += x.bytes;
            } catch (...) { std::error_code ignored; fs::remove(temporary,ignored); throw; }
        }
        x.local_trusted = std::make_unique<TrustedArchive>(destination);
        if (x.local_trusted->bytes() != x.bytes || sha256_file(x.local_trusted->pinned_path()) != x.sha)
            throw std::runtime_error("Local read-cache integrity mismatch; preserve this cache and select a new directory");
        ++local_checksum_count; check(x);
    }
    void trust(Shard & x) {
        x.trusted = std::make_unique<TrustedArchive>(dir_/x.file);
        auto started = Clock::now();
        auto digest = sha256_file(x.trusted->pinned_path());
        ++checksum_count; checksum_bytes += x.trusted->bytes(); checksum_seconds += elapsed(started);
        x.trusted->check();
        if (x.trusted->bytes() != x.bytes || digest != x.sha)
            throw std::runtime_error("Shard archive size or SHA-256 mismatch: " + x.file);
    }

    ReferenceInfo local_info(Shard const & shard) const {
        ReferenceInfo local = info_;
        local.bin_count = shard.bins;
        local.bin_insertions.assign(shard.bins,0); local.bin_distinct.assign(shard.bins,0);
        return local;
    }
    IndexMetadata metadata(Shard const & shard, ReferenceInfo const & local) const {
        auto meta = make_index_metadata(o_,local,bed_sha_);
        // Bind each format-3 archive to its canonical global range. An old
        // single archive can never pass as a shard with the same local size.
        meta.software_version += "/shard-global-start=" + std::to_string(shard.begin);
        return meta;
    }
    void checkpoint() {
        std::ostringstream payload; payload << prefix_;
        for (std::size_t id = 0; id < shards_.size(); ++id) {
            auto const & x = shards_[id];
            payload << id << '\t' << x.begin << '\t' << x.bins << '\t' << x.file << '\t'
                    << (x.sha.empty() ? "pending" : x.sha) << '\t' << x.bytes << '\t' << x.provenance << '\n';
        }
        auto text = payload.str();
        atomic_shard_text(manifest_path(),"FRAGLNCRNA-SHARDS1\n" + sha256_text(text) + "\n" + text);
    }
    void restore() {
        if (!fs::exists(manifest_path())) {
            // Refuse to adopt archives whose origin is unknown.
            for (auto const & x : shards_) if (fs::exists(dir_/x.file))
                throw std::runtime_error("Unmanaged shard archive exists without manifest; choose a new shard directory");
            checkpoint(); return;
        }
        if (fs::is_symlink(manifest_path()) || fs::file_size(manifest_path()) > 64*1024*1024)
            throw std::runtime_error("Unsafe or oversized shard manifest");
        std::ifstream input(manifest_path(),std::ios::binary);
        std::string magic, digest; std::getline(input,magic); std::getline(input,digest);
        std::string payload{std::istreambuf_iterator<char>{input},{}};
        if (magic != "FRAGLNCRNA-SHARDS1" || sha256_text(payload) != digest)
            throw std::runtime_error("Corrupt shard manifest checksum; preserve directory and rebuild in a new --ibf-shard-dir");
        if (!payload.starts_with(prefix_))
            throw std::runtime_error("Incompatible shard manifest settings/reference/mapping; choose a distinct --ibf-shard-dir for each configuration");
        std::istringstream rows(payload.substr(prefix_.size()));
        std::string line;
        for (std::size_t id = 0; id < shards_.size(); ++id) {
            if (!std::getline(rows,line)) throw std::runtime_error("Incomplete shard manifest mapping");
            auto f = split_tabs(line); auto & x = shards_[id];
            if (f.size() != 7 || f[0] != std::to_string(id) || f[1] != std::to_string(x.begin) ||
                f[2] != std::to_string(x.bins) || f[3] != x.file)
                throw std::runtime_error("Invalid shard manifest global/local mapping");
            x.sha = f[4] == "pending" ? "" : f[4]; x.bytes = parse_coordinate(f[5],"manifest archive bytes");
            x.provenance = f[6];
            if (!o_.ibf_resume && x.sha.empty()) throw std::runtime_error("Partial shard set requires --ibf-resume or a new directory");
        }
        if (std::getline(rows,line)) throw std::runtime_error("Unexpected shard manifest records");
    }
    std::unique_ptr<Ibf> build(Shard const & shard, ReferenceInfo & local, RepeatBed const & reference_bed,
                               IndexedReference const & reference) {
        auto ibf = std::make_unique<Ibf>(seqan3::bin_count{shard.bins},seqan3::bin_size{info_.bin_bits},
                                          seqan3::hash_function_count{o_.hash_functions});
        // SeqAn3 explicitly permits concurrent emplace into disjoint groups
        // of 64 local bins (interleaved storage words never share writers).
        // Each worker owns its stream, rolling encoder, BED cursor and set.
        auto groups = (shard.bins+63)/64;
        fraglncrna::parallel_records(groups,static_cast<int>(std::min(o_.actual_threads,groups)),[&](std::size_t group) {
            std::ifstream input(reference.path(),std::ios::binary);
            if (!input) throw std::runtime_error("Cannot read staged reference for shard build");
            std::array<char,65536> buffer;
            auto group_begin = shard.begin+group*64;
            auto group_end = std::min(shard.begin+shard.bins,group_begin+64);
            for (std::size_t global = group_begin; global < group_end; ++global) {
                fraglncrna::check_stop();
                auto ci = std::size_t(std::upper_bound(info_.bin_offsets.begin(),info_.bin_offsets.end(),global)-info_.bin_offsets.begin()-1);
                auto start = (global-info_.bin_offsets[ci])*o_.fragment_size;
                auto end = start + std::min(o_.fragment_size,info_.lengths[ci]-start);
                auto extended = end + std::min(o_.k-1,info_.lengths[ci]-end);
                input.clear(); input.seekg(static_cast<std::streamoff>(base_offsets_[ci]+start));
                WindowMaker maker{o_,o_.repeat == "exclude" && o_.mask_target != "query"};
                BedCursor bed{reference_bed}; bed.reset(info_.ids[ci]);
                std::unordered_set<uint64_t> distinct;
                auto local_bin = global-shard.begin;
                for (auto pos = start; pos < extended;) {
                    auto n = std::min(buffer.size(),extended-pos);
                    input.read(buffer.data(),static_cast<std::streamsize>(n));
                    if (std::size_t(input.gcount()) != n) throw std::runtime_error("Truncated staged reference in shard " + std::to_string(shard.begin));
                    for (std::size_t i = 0; i < n; ++i,++pos) {
                        auto w = maker.push(buffer[i],bed.at(pos));
                        if (!w || !w->valid) continue;
                        auto key = index_key(*w,o_); distinct.insert(key);
                        ibf->emplace(key,seqan3::bin_index{local_bin}); ++local.bin_insertions[local_bin];
                    }
                }
                local.bin_distinct[local_bin] = distinct.size();
                fraglncrna::production_log.progress("shard-build","global_bin=" + std::to_string(global) + "/" + std::to_string(info_.bin_count));
            }
        });
        return ibf;
    }
    void inspect(Shard const & shard, ReferenceInfo const & local, Ibf const & ibf) {
        std::copy(local.bin_insertions.begin(),local.bin_insertions.end(),info_.bin_insertions.begin()+shard.begin);
        std::copy(local.bin_distinct.begin(),local.bin_distinct.end(),info_.bin_distinct.begin()+shard.begin);
        auto const & raw = ibf.raw_data(); auto padded = ((shard.bins+63)/64)*64;
        for (std::size_t word = 0; word < raw.size()/64; ++word) {
            auto value = raw.data()[word]; occupied_bits += std::popcount(value);
            if (o_.mode == "calibrate") while (value) {
                auto bin = (word*64+std::countr_zero(value))%padded;
                if (bin < shard.bins) ++occupancy[shard.begin+bin];
                value &= value-1;
            }
        }
        if (o_.scoring_method == "global_ibf" &&
            (o_.legacy_global_budget || o_.genome_fpr_given || o_.ibf_auto_tune == "apply"))
            guard_bound_ += observed_ibf_union_estimate(ibf,local,o_);
    }
public:
    std::size_t built = 0, reused = 0, disk_bytes = 0, occupied_bits = 0;
    std::size_t query_loads = 0, query_bytes_read = 0;
    std::size_t local_copy_bytes = 0, local_checksum_count = 0;
    std::size_t checksum_count = 0, checksum_bytes = 0, cache_hits = 0, cache_misses = 0, cache_evictions = 0, cache_peak_bytes = 0;
    double checksum_seconds = 0, first_pass_seconds = 0, selection_seconds = 0, second_pass_seconds = 0;
    double build_seconds = 0, save_seconds = 0, load_seconds = 0, query_load_seconds = 0;
    std::vector<std::size_t> occupancy;
    ShardedIbf(Options const & o, ReferenceInfo & info, std::string bed_sha, RepeatBed const & bed,
                IndexedReference const & reference) : o_{o}, info_{info}, bed_sha_{std::move(bed_sha)},
                dir_{fs::absolute(o.ibf_shard_dir).lexically_normal()},lock_{dir_} {
        info_.bin_bits = planned_bin_bits(o.fragment_size,o.hash_functions,o.fpr);
        occupancy.resize(info.bin_count);
        base_offsets_.resize(info.ids.size());
        for (std::size_t ci = 1; ci < info.ids.size(); ++ci) base_offsets_[ci] = base_offsets_[ci-1]+info.lengths[ci-1];
        std::ostringstream spec;
        // The complete canonical compatibility JSON is readable in the manifest.
        spec << std::setprecision(17) << "{\"reference_semantic_sha256\":" << json_string(info.sequence_sha256)
             << ",\"reference_path\":" << json_string(fs::absolute(o.reference).string())
             << ",\"reference_repeat_bed_sha256\":" << json_string(bed_sha_) << ",\"k\":" << o.k
             << ",\"fragment_size\":" << o.fragment_size << ",\"strand\":" << json_string(o.strand)
             << ",\"ambiguity\":" << json_string(o.ambiguity) << ",\"repeat\":" << json_string(o.repeat)
             << ",\"mask_target\":" << json_string(o.mask_target) << ",\"p_bin_requested\":" << o.requested_per_bin_fpr
             << ",\"p_bin_effective\":" << o.fpr << ",\"hash_functions\":" << o.hash_functions
             << ",\"global_bins\":" << info.bin_count << ",\"bin_bits\":" << info.bin_bits
             << ",\"max_bins_per_ibf\":" << o.max_bins_per_ibf << ",\"compiler\":" << json_string(__VERSION__)
             << ",\"seqan\":\"3.4.1\",\"index_format\":3,\"size_t_bytes\":" << sizeof(std::size_t)
             << ",\"byte_order\":" << json_string(std::endian::native == std::endian::little ? "little" : "big") << "}\n";
        spec << "global_bin\tshard_id\tlocal_bin\tcontig\tstart\tend\textended_end\n";
        for (std::size_t ci = 0; ci < info.ids.size(); ++ci)
            for (auto start = std::size_t{0}; start < info.lengths[ci];) {
                auto global = info.bin_offsets[ci]+start/o.fragment_size;
                auto end = start+std::min(o.fragment_size,info.lengths[ci]-start);
                spec << global << '\t' << global/o.max_bins_per_ibf << '\t' << global%o.max_bins_per_ibf
                     << '\t' << info.ids[ci] << '\t' << start << '\t' << end << '\t'
                     << end+std::min(o.k-1,info.lengths[ci]-end) << '\n'; start = end;
            }
        spec << "shard_id\tglobal_begin\treal_bins\trelative_path\tsha256\tbytes\tbuild_run\n";
        prefix_ = spec.str();
        for (std::size_t begin = 0; begin < info.bin_count;) {
            Shard x; x.begin = begin; x.bins = std::min(o.max_bins_per_ibf,info.bin_count-begin);
            x.file = "shard-" + std::to_string(shards_.size()) + ".ibf"; shards_.push_back(std::move(x));
            begin += shards_.back().bins;
        }
        // Bound cache by configured bytes, 1/32 physical RAM, and explicit
        // resident-index budget. A transient uncached load also fits the cap.
        auto physical = double(fraglncrna::memory_limit_bytes());
        cache_budget_ = static_cast<std::size_t>(std::min(double(o.ibf_shard_cache_mb)*1048576,physical/32));
        if (o.ibf_max_index_memory_gb > 0)
            cache_budget_ = std::min(cache_budget_,static_cast<std::size_t>(o.ibf_max_index_memory_gb*1e9));
        fraglncrna::production_log.event("shard-cache","stage_completed","budget_bytes=" + std::to_string(cache_budget_) +
            " policy=deterministic_LRU integrity=SHA256_once_pinned_descriptor_and_change_stamps");
        if (dir_.string().starts_with("/mnt/"))
            fraglncrna::production_log.event("shard-cache","warning","Mounted-volume index IO may be slow under WSL; use an explicitly staged verified index copy on native Linux storage for comparison. Canonical index path is unchanged.",3);
        restore();
        if (!o.ibf_local_cache_dir.empty()) {
            if (fs::is_symlink(o.ibf_local_cache_dir)) throw std::runtime_error("Symlink local read-cache directory refused");
            local_root_ = fs::absolute(o.ibf_local_cache_dir);
            local_lock_ = std::make_unique<ShardDirectoryLock>(local_root_);
            local_dir_ = local_root_/sha256_text(prefix_);
            if (fs::is_symlink(local_dir_)) throw std::runtime_error("Symlink local cache compatibility directory refused");
            fs::create_directories(local_dir_);
            auto projected = planned_index_bytes(std::min(o.max_bins_per_ibf,info.bin_count),info.bin_bits)*shards_.size();
            if (projected > o.ibf_local_cache_max_mb*1048576)
                throw std::runtime_error("Explicit local read-cache archive set exceeds --ibf-local-cache-max-mb");
            auto existing = local_used_bytes();
            if (existing > o.ibf_local_cache_max_mb*1048576)
                throw std::runtime_error("Existing local cache exceeds byte budget; choose a new cache directory");
            fraglncrna::production_log.event("shard-cache","stage_started","explicit_local_read_cache=" + local_dir_.string() +
                " projected_archive_bytes=" + std::to_string(projected) + " canonical_index=" + dir_.string());
        }
        for (std::size_t id = 0; id < shards_.size(); ++id) {
            fraglncrna::check_stop(); auto & x = shards_[id]; auto path = dir_/x.file;
            if (fs::is_symlink(path)) throw std::runtime_error("Refusing symlink shard archive " + path.string());
            auto local = local_info(x); auto expected = metadata(x,local);
            auto load_before = load_seconds, build_before = build_seconds, save_before = save_seconds;
            std::string outcome = x.sha.empty() ? "built" : "rebuilt";
            std::unique_ptr<Ibf> ibf;
            make_room(x);
            if (!x.sha.empty() && o.ibf_shard_reuse) {
                auto started = Clock::now();
                try {
                    trust(x);
                    ibf = load_index_file(read_path(x),expected,local,o.ibf_max_index_memory_gb,true);
                    x.trusted->check(); ++reused; outcome = "reused";
                } catch (std::exception const & e) {
                    if (!o.ibf_resume) throw;
                    fraglncrna::production_log.event("shard","warning","invalidated shard=" + std::to_string(id) + " reason=" + e.what(),3);
                }
                load_seconds += elapsed(started);
            }
            if (!ibf) {
                if (!o.ibf_resume && !x.sha.empty() && o.ibf_shard_reuse)
                    throw std::runtime_error("Invalid shard requires --ibf-resume");
                fraglncrna::production_log.event("shard","stage_started","build shard=" + std::to_string(id) +
                    " real_bins=" + std::to_string(x.bins) + " insertion_threads=" + std::to_string(std::min(o_.actual_threads,(x.bins+63)/64)));
                auto started = Clock::now(); ibf = build(x,local,bed,reference); build_seconds += elapsed(started);
                x.owned = !fs::exists(path); x.provenance = fraglncrna::production_log.run_id();
                started = Clock::now();
                save_index_file(path,metadata(x,local),*ibf);
                // Verify complete serialized archive before advertising it in a checkpoint.
                if (o.ibf_shard_storage == "disk") ibf.reset();
                int fd = open(path.c_str(),O_RDONLY | O_NOFOLLOW);
                if (fd < 0) throw std::runtime_error("Cannot flush shard archive");
                auto flushed = fsync(fd); close(fd);
                if (flushed) throw std::runtime_error("Cannot flush shard archive");
                sync_shard_directory(dir_);
                x.trusted = std::make_unique<TrustedArchive>(path);
                auto hash_start = Clock::now();
                x.sha = sha256_file(x.trusted->pinned_path()); x.bytes = x.trusted->bytes();
                ++checksum_count; checksum_bytes += x.bytes; checksum_seconds += elapsed(hash_start);
                x.trusted->check();
                if (o.ibf_shard_storage == "disk")
                    // A newly saved file has no trusted manifest digest yet.
                    // Check its independent embedded payload checksum before
                    // advertising the freshly computed whole-file digest.
                    ibf = load_index_file(x.trusted->pinned_path(),expected,local,o.ibf_max_index_memory_gb);
                else {
                    // Validate durable archive checksums/metadata while retaining
                    // the live built object. Querying never deserializes it again.
                    auto pinned = x.trusted->pinned_path();
                    std::ifstream archive(pinned,std::ios::binary);
                    auto stored = validate_index_archive(pinned,metadata(x,local),archive);
                    if (stored.bin_bits != ibf->bin_size() || stored.bin_count != ibf->bin_count())
                        throw std::runtime_error("Saved shard dimensions differ from live built IBF");
                }
                x.trusted->check(); save_seconds += elapsed(started);
                ++built; checkpoint();
            }
            if (!x.local_trusted) stage_local(x);
            auto occupied_before = occupied_bits;
            inspect(x,local,*ibf); disk_bytes += x.bytes;
            fraglncrna::production_log.event("shard","stage_completed","shard=" + std::to_string(id) +
                " status=" + outcome + " current_run_owned=" + std::to_string(x.owned) +
                " bytes=" + std::to_string(x.bytes) + " sha256=" + x.sha +
                " inserted_windows=" + std::to_string(std::accumulate(local.bin_insertions.begin(),local.bin_insertions.end(),std::size_t{0})) +
                " distinct_keys=" + std::to_string(std::accumulate(local.bin_distinct.begin(),local.bin_distinct.end(),std::size_t{0})) +
                " occupied_bits=" + std::to_string(occupied_bits-occupied_before) +
                " build_seconds=" + std::to_string(build_seconds-build_before) +
                " load_seconds=" + std::to_string(load_seconds-load_before) +
                " save_and_validate_seconds=" + std::to_string(save_seconds-save_before) +
                " bin_bits=" + std::to_string(info.bin_bits) + " padded_bins=" + std::to_string(((x.bins+63)/64)*64) +
                " peak_RSS=" + std::to_string(peak_rss_bytes()));
            retain(x,std::move(ibf));
            // Disk-mode local IBF is destroyed before next allocation.
        }
    }
    fs::path manifest_path() const { return dir_/"manifest.tsv"; }
    std::size_t size() const { return shards_.size(); }
    std::size_t max_resident() const { return resident_peak_; }
    std::size_t cache_budget_bytes() const { return cache_budget_; }
    void check_integrity() const { for (auto const & x : shards_) check(x); }
    double occupancy_union_estimate() const { return std::min(1.0,guard_bound_); }
    std::size_t begin(std::size_t id) const { return shards_.at(id).begin; }
    std::size_t bins(std::size_t id) const { return shards_.at(id).bins; }
    template <typename Fn> void visit(std::size_t id, Fn fn) {
        auto & x = shards_.at(id); check(x);
        x.last_use = ++cache_clock_;
        if (x.memory) { ++cache_hits; fn(*x.memory); check(x); return; }
        ++cache_misses;
        make_room(x);
        auto local = local_info(x); auto started = Clock::now();
        auto ibf = load_index_file(read_path(x),metadata(x,local),local,o_.ibf_max_index_memory_gb,true);
        check(x);
        query_load_seconds += elapsed(started); ++query_loads; query_bytes_read += x.bytes;
        std::size_t count = 1; for (auto const & other : shards_) count += bool(other.memory);
        resident_peak_ = std::max(resident_peak_,count);
        fn(*ibf); check(x); retain(x,std::move(ibf));
    }
    void cleanup_success() {
        if (o_.ibf_shard_retention != "cleanup_on_success") return;
        for (auto & x : shards_) if (x.owned) {
            auto path = dir_/x.file;
            if (!fs::is_symlink(path) && fs::exists(path) && sha256_file(path) == x.sha) {
                fs::remove(path);
                fraglncrna::production_log.event("cleanup","stage_completed","removed owned shard=" + path.string());
            }
        }
        // Keep the authoritative manifest and ownership provenance. Missing
        // owned archives are explicit and safely rebuildable on the next run.
    }
};
