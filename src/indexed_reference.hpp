#pragma once

// Reusable random-access cache of original FASTA bases, including case/Ns.
// Metadata remains in ReferenceInfo; contigs are never searched as a joined
// sequence. Cache content is checked against the freshly inspected reference.
class IndexedReference {
    fs::path path_, temporary_;
    bool ephemeral_ = false, reused_ = false;
    std::ofstream output_;
    std::string buffer_;
    fraglncrna::Sha256 digest_;
    std::size_t bases_ = 0;
    struct Page { std::string data; std::list<std::size_t>::iterator recency; };
    mutable std::mutex page_mutex_;
    mutable std::unordered_map<std::size_t,Page> pages_;
    mutable std::list<std::size_t> recency_;
    std::size_t page_limit_ = 0;
    static constexpr std::size_t page_size_ = 16384;
    mutable std::size_t page_hits_ = 0, page_misses_ = 0, read_bytes_ = 0;
    mutable double read_seconds_ = 0, lock_wait_seconds_ = 0;
public:
    explicit IndexedReference(fs::path requested, std::size_t page_cache_mb = 64) {
        auto physical = double(fraglncrna::memory_limit_bytes());
        page_limit_ = static_cast<std::size_t>(std::min(double(page_cache_mb)*1048576,physical/64))/page_size_;
        if (requested.empty()) {
            std::string pattern = (fs::temp_directory_path() / "fraglncrna-reference-XXXXXX").string();
            std::vector<char> name(pattern.begin(), pattern.end()); name.push_back('\0');
            int fd = mkstemp(name.data());
            if (fd < 0) throw std::runtime_error("Cannot create indexed reference cache in temporary directory");
            close(fd); path_ = name.data(); ephemeral_ = true;
        } else path_ = std::move(requested);
        reused_ = !ephemeral_ && fs::exists(path_);
        if (!reused_) {
            fs::create_directories(path_.parent_path().empty() ? "." : path_.parent_path());
            temporary_ = unique_staging_path(path_);
            output_.open(temporary_, std::ios::binary);
            if (!output_) throw std::runtime_error("Cannot create indexed reference cache " + temporary_.string());
            buffer_.reserve(65536);
        }
    }
    ~IndexedReference() {
        std::error_code ignored;
        if (!temporary_.empty()) fs::remove(temporary_, ignored);
        if (ephemeral_) fs::remove(path_, ignored);
    }
    void base(char c) {
        digest_.update(static_cast<unsigned char>(c)); ++bases_;
        if (!reused_) {
            buffer_ += c;
            if (buffer_.size() == 65536) flush();
        }
    }
    void flush() {
        output_.write(buffer_.data(), static_cast<std::streamsize>(buffer_.size())); buffer_.clear();
        if (!output_) throw std::runtime_error("Indexed reference cache write failed; check free disk space");
    }
    void finish() {
        if (reused_) {
            if (fs::file_size(path_) != bases_ || sha256_file(path_) != digest_.hex_digest())
                throw std::runtime_error("Indexed reference cache checksum mismatch: choose a new --reference-cache path");
        } else {
            flush(); output_.close();
            if (!output_) throw std::runtime_error("Cannot finalize indexed reference cache");
            fs::rename(temporary_, path_); temporary_.clear();
        }
        fraglncrna::production_log.event("reference-cache", "stage_completed",
            "strategy=normalized_original_bases path=" + path_.string() + " reused=" + std::to_string(reused_) +
            " preprocessing_bases=" + std::to_string(bases_) + " cache_sha256=" + digest_.hex_digest());
    }
    // The lock serializes page fill/eviction and copying, so borrowed page
    // pointers never outlive an eviction. Decoding proceeds outside the lock.
    void read(std::ifstream & input, std::size_t offset, char * target, std::size_t length) const {
        if (offset > bases_ || length > bases_-offset) throw std::runtime_error("Reference page request outside validated bounds");
        auto waiting = Clock::now();
        std::lock_guard lock(page_mutex_);
        lock_wait_seconds_ += std::chrono::duration<double>(Clock::now()-waiting).count();
        while (length) {
            auto page = offset/page_size_, within = offset%page_size_;
            auto available = std::min(page_size_,bases_-page*page_size_);
            auto count = std::min(length,available-within);
            if (!count) throw std::runtime_error("Reference page outside validated bounds");
            if (!page_limit_) {
                auto reading = Clock::now();
                input.clear(); input.seekg(static_cast<std::streamoff>(offset));
                input.read(target,static_cast<std::streamsize>(count));
                if (std::size_t(input.gcount()) != count) throw std::runtime_error("Truncated reference page");
                read_seconds_ += std::chrono::duration<double>(Clock::now()-reading).count();
                ++page_misses_; read_bytes_ += count;
            } else {
                auto it = pages_.find(page);
                if (it == pages_.end()) {
                    if (pages_.size() == page_limit_) { pages_.erase(recency_.back()); recency_.pop_back(); }
                    std::string data(available,'\0');
                    auto reading = Clock::now();
                    input.clear(); input.seekg(static_cast<std::streamoff>(page*page_size_));
                    input.read(data.data(),static_cast<std::streamsize>(data.size()));
                    if (std::size_t(input.gcount()) != data.size()) throw std::runtime_error("Truncated reference page");
                    read_seconds_ += std::chrono::duration<double>(Clock::now()-reading).count();
                    ++page_misses_; read_bytes_ += data.size(); recency_.push_front(page);
                    it = pages_.emplace(page,Page{std::move(data),recency_.begin()}).first;
                } else { ++page_hits_; recency_.splice(recency_.begin(),recency_,it->second.recency); }
                std::copy_n(it->second.data.data()+within,count,target);
            }
            offset += count; target += count; length -= count;
        }
    }
    double read_seconds() const { return read_seconds_; }
    double lock_wait_seconds() const { return lock_wait_seconds_; }
    std::size_t page_hits() const { return page_hits_; }
    std::size_t page_misses() const { return page_misses_; }
    std::size_t read_bytes() const { return read_bytes_; }
    std::size_t page_cache_bytes() const { return page_limit_*page_size_; }
    fs::path const & path() const { return path_; }
    bool reused() const { return reused_; }
};
