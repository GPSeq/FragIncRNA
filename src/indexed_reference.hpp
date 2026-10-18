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
public:
    explicit IndexedReference(fs::path requested) {
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
    fs::path const & path() const { return path_; }
    bool reused() const { return reused_; }
};
