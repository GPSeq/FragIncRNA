#pragma once

// Optional reports stream in original query/window order. No window-by-bin
// matrix is retained. Both gzip and plain outputs are staged until success.
class ReportStream {
    fs::path final_, staged_;
    std::ofstream plain_;
    gzFile gzip_ = nullptr;
public:
    explicit ReportStream(fs::path const & path) : final_{path} {
        if (path.empty()) return;
        auto parent = path.parent_path().empty() ? fs::path{"."} : path.parent_path();
        fs::create_directories(parent);
        if (fs::space(parent).available < 1048576) throw std::runtime_error("Insufficient disk space for window report: " + path.string());
        staged_ = unique_staging_path(path);
        if (path.extension() == ".gz") {
            gzip_ = gzopen(staged_.c_str(),"wb");
            if (!gzip_) throw std::runtime_error("Cannot open compressed report " + path.string());
        } else {
            plain_.open(staged_,std::ios::binary);
            if (!plain_) throw std::runtime_error("Cannot open report " + path.string());
        }
    }
    ~ReportStream() { if (gzip_) gzclose(gzip_); }
    void write(std::string const & text) {
        if (gzip_) {
            if (gzwrite(gzip_,text.data(),static_cast<unsigned>(text.size())) != static_cast<int>(text.size()))
                throw std::runtime_error("Window report gzip write failed; check disk space");
        } else if (plain_.is_open()) {
            plain_ << text; if (!plain_) throw std::runtime_error("Window report write failed; check disk space");
        }
    }
    void finish() {
        if (final_.empty()) return;
        if (gzip_) { auto code = gzclose(gzip_); gzip_ = nullptr; if (code != Z_OK) throw std::runtime_error("Cannot finalize gzip report"); }
        if (plain_.is_open()) { plain_.close(); if (!plain_) throw std::runtime_error("Cannot finalize report"); }
        fs::rename(staged_,final_); staged_.clear();
    }
};

std::string fasta_safe_id(std::string const & text) {
    std::ostringstream out;
    for (unsigned char c : text) {
        if (std::isalnum(c) || c == '_' || c == '-' || c == '.') out << c;
        else out << '%' << std::uppercase << std::hex << std::setw(2) << std::setfill('0') << unsigned(c);
    }
    return out.str();
}

class WindowReport {
    Options const & o_;
    ReportStream table_, fasta_;
    std::unordered_set<std::string> selected_, remaining_;
    bool filter_ = false;
    std::size_t rows_ = 0;
public:
    explicit WindowReport(Options const & o) : o_{o}, table_{o.window_output}, fasta_{o.kmer_fasta_output} {
        if (!o.window_query_ids.empty()) {
            filter_ = true; std::ifstream in(o.window_query_ids); std::string id;
            if (!in) throw std::runtime_error("Cannot open --window-query-ids");
            while (std::getline(in,id)) { if (!id.empty() && id.back() == '\r') id.pop_back(); if (!id.empty()) selected_.insert(id); }
            remaining_ = selected_;
        }
        if (!o.window_output.empty()) table_.write("schema_version\treference_id\tk\tquery_index\tquery_id\tk_index\tquery_start\tquery_end\tvalid\tcanonical_kmer_id\tibf_any_bin_positive\tibf_candidate_positive\texact_verified_in_candidates\tmatching_candidate_regions_count" +
            std::string(o.window_include_sequence ? "\tkmer_sequence" : "") +
            "\tscoring_method\tstrand_mode\trepeat_policy\tmask_target\tverification_scope\n");
        if (!o.window_output.empty() || !o.kmer_fasta_output.empty())
            std::cerr << "Window reports enabled: scope=" << o.window_scope << " row_cap=" << o.max_window_rows <<
                " (all exports every sliding position, including invalid windows; limit failure leaves staged files). "
                << "Rough TSV size at the cap: " << double(o.max_window_rows)*256 << " bytes; actual size depends on IDs/sequences.\n";
    }
    void write(std::string const & reference_id, std::vector<Query> const & batch, ExactResult const & exact) {
        if (o_.window_output.empty() && o_.kmer_fasta_output.empty()) return;
        bool queried = o_.mode != "exact" && o_.mode != "export-kmers";
        bool verified = o_.mode != "export-kmers" && o_.mode != "exact" && o_.fdr_control_per_bin && o_.per_bin_exact_verify;
        for (std::size_t qi = 0; qi < batch.size(); ++qi) {
            auto const & q = batch[qi]; remaining_.erase(q.id);
            if (filter_ && !selected_.contains(q.id)) continue;
            for (std::size_t wi = 0; wi < q.windows.size(); ++wi) {
                auto const & w = q.windows[wi];
                bool candidate = o_.fdr_control_per_bin && q.candidate_support.size() > wi && q.candidate_support[wi];
                bool matched = verified && exact.candidate_window_masks.size() > qi &&
                    (exact.candidate_window_masks[qi][wi/64] & (uint64_t{1} << (wi%64)));
                if (o_.window_scope == "verified" && !(verified && w.valid && matched)) continue;
                if (o_.window_scope == "candidate" && !(queried && w.valid && candidate)) continue;
                if (o_.window_scope == "failed" && !(verified && w.valid && !matched)) continue;
                if (++rows_ > o_.max_window_rows) throw std::runtime_error("Window report exceeds --max-window-rows query=" + q.id + " k_index=" + std::to_string(wi));
                std::ostringstream canonical;
                if (w.valid) canonical << "dna2bit-v1:" << o_.k << ':' << std::hex << std::min(w.forward,w.reverse);
                else canonical << "NA";
                auto sequence = w.valid && !q.report_sequence.empty() ? q.report_sequence.substr(wi,o_.k) : "NA";
                std::ostringstream row;
                row << "1\t" << reference_id << '\t' << o_.k << '\t' << q.index << '\t' << q.id << '\t' << wi << '\t' << wi << '\t' << wi+o_.k
                    << '\t' << w.valid << '\t' << canonical.str() << '\t'
                    << (queried && w.valid ? (q.bin_hits[wi] ? "1" : "0") : "NA") << '\t'
                    << (queried && o_.fdr_control_per_bin && w.valid ? (candidate ? "1" : "0") : "NA") << '\t'
                    << (verified && w.valid ? (matched ? "1" : "0") : "NA") << '\t';
                if (verified && w.valid) row << (exact.candidate_region_matches.size() > qi ? exact.candidate_region_matches[qi][wi] : 0);
                else row << "NA";
                if (o_.window_include_sequence) row << '\t' << sequence;
                row << '\t' << (o_.mode == "export-kmers" ? "query_only" : o_.scoring_method) << '\t' << o_.strand << '\t' << o_.repeat << '\t' << o_.mask_target
                    << '\t' << (verified ? "candidate_intervals" : "unavailable") << '\n';
                table_.write(row.str());
                if (w.valid && !o_.kmer_fasta_output.empty())
                    fasta_.write(">ref=" + fasta_safe_id(reference_id) + "|qi=" + std::to_string(q.index) + "|q=" + fasta_safe_id(q.id) +
                        "|k=" + std::to_string(o_.k) + "|wi=" + std::to_string(wi) + "|start=" + std::to_string(wi) + "|end=" + std::to_string(wi+o_.k) + "\n" + sequence + "\n");
            }
        }
    }
    std::size_t rows() const { return rows_; }
    void finish() {
        if (!remaining_.empty()) throw std::runtime_error("Window query ID absent from input FASTA: " + *remaining_.begin());
        table_.finish(); fasta_.finish();
    }
};
