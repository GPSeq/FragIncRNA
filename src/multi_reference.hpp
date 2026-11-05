#pragma once

std::string lower_ascii(std::string text) {
    for (char & c : text) if (c >= 'A' && c <= 'Z') c += 'a'-'A';
    return text;
}
std::string directory_reference_id(fs::path const & file) {
    auto name = file.filename().string(), lower = lower_ascii(name);
    for (auto const & suffix : {".fasta.gz",".fna.gz",".fa.gz",".fasta",".fna",".fa"})
        if (lower.ends_with(suffix)) {
            name.resize(name.size()-std::char_traits<char>::length(suffix));
            for (char & c : name) if (!((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
                (c >= '0' && c <= '9') || c == '_' || c == '-' || c == '.')) c = '_';
            if (name.empty() || name == "." || name == ".." || name.front() == '.')
                throw std::runtime_error("Unsafe/empty reference ID: " + file.string());
            return name;
        }
    return {};
}
std::string replace_token(std::string text, std::string const & token, std::string const & value) {
    for (std::size_t at = 0; (at = text.find(token,at)) != std::string::npos; at += value.size()) text.replace(at,token.size(),value);
    return text;
}

int run_reference_directory(Options const & o) {
    if (!o.save_index.empty() || !o.load_index.empty() || o.ibf_build_only || o.ibf_dry_run)
        throw std::runtime_error("Directory mode does not accept single save/load archives, build-only or dry-run; use --reference FILE");
    if (o.ibf_shard_retention != "keep") throw std::runtime_error("Directory resume requires ibf_shard_retention=keep");
    if (!fs::is_directory(o.reference_dir) || fs::is_symlink(o.reference_dir))
        throw std::runtime_error("--reference-dir must be a real directory (no symlink)");
    auto out_root = fs::absolute(o.output_dir).lexically_normal();
    if (fs::is_symlink(out_root)) throw std::runtime_error("Symlink output root refused");
    fs::create_directories(out_root);
    ShardDirectoryLock run_lock{out_root};
    auto index_root = o.ibf_shard_dir.empty() ? out_root/".indexes" : fs::absolute(o.ibf_shard_dir).lexically_normal();
    struct Entry {
        fs::path path, output, index;
        std::string id, hash, status = "PENDING", validation = "unprocessed", error;
        double seconds = 0; std::size_t queries = 0;
    };
    std::vector<Entry> entries;
    auto discover = [&](auto const & item) {
        if (item.is_symlink() || !item.is_regular_file()) return;
        auto id = directory_reference_id(item.path()); if (id.empty()) return;
        auto path = fs::absolute(item.path()).lexically_normal();
        Entry entry; entry.path = path; entry.id = id;
        entry.output = out_root/id/("k"+std::to_string(o.k)); entry.index = index_root/id/("k"+std::to_string(o.k));
        entries.push_back(std::move(entry));
    };
    if (o.reference_recursive) for (auto const & item : fs::recursive_directory_iterator(o.reference_dir)) discover(item);
    else for (auto const & item : fs::directory_iterator(o.reference_dir)) discover(item);
    std::sort(entries.begin(),entries.end(),[](auto const & a,auto const & b) { return a.path.string() < b.path.string(); });
    if (entries.empty()) throw std::runtime_error("No recognized FASTA files in --reference-dir");
    std::unordered_set<std::string> ids;
    for (auto const & e : entries) if (!ids.insert(lower_ascii(e.id)).second)
        throw std::runtime_error("Case-insensitive/sanitized reference ID collision: " + e.id);
    std::vector<std::unique_ptr<TrustedArchive>> shared_inputs;
    for (auto const & input : {o.queries,o.config_path,o.threshold_table,o.reference_repeat_bed,o.query_repeat_bed,o.exon_blocks,o.expected_loci_bed,o.window_query_ids})
        if (!input.empty()) shared_inputs.push_back(std::make_unique<TrustedArchive>(input));
    auto check_inputs = [&] { for (auto const & input : shared_inputs) input->check(); };
    check_inputs();
    // Conservative invocation fingerprint: changing even a resource option
    // invalidates resume rather than risk accepting incompatible results.
    std::ostringstream settings;
    for (auto const & arg : o.expanded_arguments) settings << arg.size() << ':' << arg << '\n';
    settings << sha256_file(o.config_path) << sha256_file(o.queries) << sha256_file("/proc/self/exe");
    for (auto const & input : {o.threshold_table,o.reference_repeat_bed,o.query_repeat_bed,o.exon_blocks,o.expected_loci_bed,o.window_query_ids})
        settings << sha256_file(input);
    auto settings_sha = sha256_text(settings.str());
    for (auto & e : entries) e.hash = sha256_file(e.path);
    auto manifest = out_root/"run_manifest.tsv";
    if (fs::exists(manifest)) {
        if (fs::is_symlink(manifest) || !fs::is_regular_file(manifest)) throw std::runtime_error("Unsafe run manifest");
        std::ifstream old(manifest); std::string header; std::getline(old,header);
        if (!header.starts_with("schema_version\treference_id\treference_path\treference_file_sha256\tk\tsettings_sha256\t"))
            throw std::runtime_error("Unrecognized run manifest; use a new --output-dir");
    }
    auto publish_manifest = [&] {
        std::ostringstream text;
        text << "schema_version\treference_id\treference_path\treference_file_sha256\tk\tsettings_sha256\tindex_directory\toutput_directory\tstatus\truntime_seconds\tquery_count\tvalidation\terror\n";
        for (auto const & e : entries)
            text << "1\t" << e.id << '\t' << e.path.string() << '\t' << e.hash << '\t' << o.k << '\t' << settings_sha << '\t'
                << e.index.string() << '\t' << e.output.string() << '\t' << e.status << '\t' << e.seconds << '\t' << e.queries << '\t' << e.validation << '\t' << e.error << '\n';
        atomic_shard_text(manifest,text.str());
    };
    // An assembly directory contains a sealed ledger of its files. A crash
    // after atomic directory publication but before manifest update is safe:
    // validate the ledger and recover the completed state on the next run.
    auto validate_complete = [&](Entry & e) {
        auto seal = e.output/"COMPLETED.sha256.tsv";
        if (!fs::is_regular_file(seal) || fs::is_symlink(seal)) throw std::runtime_error("Existing output has no completion ledger: " + e.output.string());
        std::ifstream in(seal); std::string header;
        std::getline(in,header);
        auto seal_header = split_tabs(header);
        if (seal_header.size() != 5 || seal_header[0] != "FRAGLNCRNA-COMPLETED1" || seal_header[1] != settings_sha || seal_header[2] != e.hash)
            throw std::runtime_error("Completed output settings/reference changed; use a new --output-dir: " + e.id);
        e.seconds = std::stod(seal_header[3]); e.queries = std::stoull(seal_header[4]);
        std::size_t count = 0; std::string line;
        while (std::getline(in,line)) {
            auto fields = split_tabs(line);
            if (fields.size() != 2 || fs::path(fields[1]).filename() != fields[1] || fields[1] == "COMPLETED.sha256.tsv")
                throw std::runtime_error("Invalid completion ledger");
            auto file = e.output/fields[1];
            if (fs::is_symlink(file) || !fs::is_regular_file(file) || sha256_file(file) != fields[0])
                throw std::runtime_error("Completed output integrity mismatch: " + file.string());
            ++count;
        }
        if (count < 2 || !in.eof()) throw std::runtime_error("Incomplete completion ledger");
        e.status = "COMPLETED"; e.validation = "completion_ledger_verified";
    };
    for (auto & e : entries) {
        if (e.path.string().find_first_of("\t\r\n") != std::string::npos || e.output.string().find_first_of("\t\r\n") != std::string::npos)
            throw std::runtime_error("Directory run paths must not contain tabs/newlines");
        if (fs::is_symlink(e.output) || fs::is_symlink(e.output.parent_path())) throw std::runtime_error("Symlink output directory refused");
        if (fs::exists(e.output)) validate_complete(e);
    }
    publish_manifest();
    bool failed = false;
    for (std::size_t ri = 0; ri < entries.size(); ++ri) {
        auto & e = entries[ri];
        if (e.status == "COMPLETED") { std::cerr << "Genome " << ri+1 << '/' << entries.size() << " resume verified " << e.id << '\n'; continue; }
        std::cerr << "Genome " << ri+1 << '/' << entries.size() << " start " << e.id << '\n';
        fs::create_directories(e.output.parent_path());
        std::string pattern = (e.output.parent_path()/".k-staging-XXXXXX").string();
        std::vector<char> name(pattern.begin(),pattern.end()); name.push_back('\0');
        auto made = mkdtemp(name.data()); if (!made) throw std::runtime_error("Cannot create genome staging directory");
        fs::path stage{made}; auto started = Clock::now();
        e.status = "RUNNING"; e.validation = "staged"; publish_manifest();
        try {
            check_inputs();
            std::vector<std::string> args;
            std::unordered_set<std::string> scalar_drop{"--reference-dir","--output-dir","--output-name-template","--ibf-shard-dir","--reference-cache","--reference-id"};
            std::unordered_set<std::string> flag_drop{"--reference-recursive","--no-reference-recursive","--continue-on-reference-error","--no-continue-on-reference-error"};
            std::unordered_set<std::string> sidecars{"--metrics","--log-file","--ibf-tuning-report","--positions","--calibration-bins","--fragments","--repeat-classes","--reference-hit-classes","--discover-loci","--expected-loci-output","--per-bin-output","--window-output","--kmer-fasta-output"};
            for (std::size_t i = 0; i < o.expanded_arguments.size(); ++i) {
                auto const & arg = o.expanded_arguments[i];
                if (scalar_drop.contains(arg)) { ++i; continue; }
                if (flag_drop.contains(arg)) continue;
                args.push_back(arg);
                if (sidecars.contains(arg)) {
                    auto value = fs::path{o.expanded_arguments.at(++i)};
                    if (!value.empty() && (value.filename() != value || value == "." || value == ".." || lower_ascii(value.string()) == "completed.sha256.tsv"))
                        throw std::runtime_error("Directory sidecar options require safe non-reserved basenames: " + arg);
                    args.push_back(value.empty() ? "" : (stage/value).string());
                }
            }
            auto result_name = replace_token(replace_token(o.output_name_template,"{reference_id}",e.id),"{k}",std::to_string(o.k));
            if (lower_ascii(result_name) == "completed.sha256.tsv") throw std::runtime_error("Reserved output filename");
            if (result_name.empty() || fs::path(result_name).filename() != result_name || result_name == "." || result_name == ".." || result_name.find_first_of("{}\t\r\n") != std::string::npos)
                throw std::runtime_error("--output-name-template must expand to a safe basename");
            auto add = [&](std::string key,fs::path const & value) { args.push_back(std::move(key)); args.push_back(value.string()); };
            add("--reference",e.path); add("--reference-id",e.id); add("--output",stage/result_name);
            if (o.ibf_layout != "single" && o.mode != "exact") add("--ibf-shard-dir",e.index);
            if (!o.reference_cache.empty()) add("--reference-cache",fs::absolute(o.reference_cache)/e.id/("k"+std::to_string(o.k)+".bases"));
            // All normal parse/alias/input/output checks apply to each genome.
            std::vector<char *> argv; for (auto & arg : args) argv.push_back(arg.data());
            auto child = parse_options(static_cast<int>(argv.size()),argv.data()); child.config_path = o.config_path;
            if (run(child)) throw std::runtime_error("Genome returned nonzero status");
            fraglncrna::production_log.close();
            // Never publish if original reference or shared inputs changed.
            check_inputs();
            if (sha256_file(e.path) != e.hash) throw std::runtime_error("Reference changed during directory run");
            e.seconds = elapsed(started); e.status = "COMPLETED"; e.validation = "completed_atomic_directory";
            // Cheap completed query count from the main transcript TSV.
            std::ifstream rows(stage/result_name); std::string line; std::getline(rows,line);
            while (std::getline(rows,line)) ++e.queries;
            // Stored provenance names the final publication directory. The
            // only mutable copies are this run's still-private text sidecars.
            for (auto const & file : {child.metrics,child.log_file,fs::path(child.log_file.string()+".jsonl"),child.ibf_tuning_report}) {
                if (!fs::exists(file)) continue;
                std::ifstream in(file,std::ios::binary); std::string text{std::istreambuf_iterator<char>{in},{}}; in.close();
                if (file.extension() == ".json" || file.extension() == ".jsonl")
                    text = replace_token(text,fraglncrna::log_escape(stage.string()),fraglncrna::log_escape(e.output.string()));
                else text = replace_token(text,stage.string(),e.output.string());
                std::ofstream out(file,std::ios::binary | std::ios::trunc); out << text; out.close();
                if (!out) throw std::runtime_error("Cannot finalize directory provenance");
            }
            std::vector<fs::path> files;
            for (auto const & f : fs::directory_iterator(stage)) if (f.is_regular_file() && !f.is_symlink()) files.push_back(f.path());
            std::sort(files.begin(),files.end());
            std::ostringstream ledger; ledger << "FRAGLNCRNA-COMPLETED1\t" << settings_sha << '\t' << e.hash << '\t' << e.seconds << '\t' << e.queries << '\n';
            for (auto const & file : files) {
                ledger << sha256_file(file) << '\t' << file.filename().string() << '\n';
                int descriptor = open(file.c_str(),O_RDONLY | O_NOFOLLOW | O_CLOEXEC);
                if (descriptor < 0) throw std::runtime_error("Cannot open completed output for sync");
                auto result = fsync(descriptor); close(descriptor);
                if (result) throw std::runtime_error("Cannot sync completed output");
            }
            atomic_shard_text(stage/"COMPLETED.sha256.tsv",ledger.str());
            if (fs::exists(e.output)) throw std::runtime_error("Output appeared during run; refusing replacement");
            fs::rename(stage,e.output); sync_shard_directory(e.output.parent_path()); publish_manifest();
        } catch (std::exception const & error) {
            fraglncrna::production_log.close(); e.seconds = elapsed(started); e.status = "FAILED"; e.validation = "partial_staging_not_published";
            e.error = error.what(); for (auto & c : e.error) if (c == '\t' || c == '\r' || c == '\n') c = ' ';
            publish_manifest(); std::cerr << "Genome failed " << e.id << ": " << e.error << " staging=" << stage.string() << '\n';
            failed = true;
            if (!o.continue_on_reference_error || fraglncrna::requested_signal.load()) throw;
        }
    }
    return failed ? 1 : 0;
}
