#pragma once

// Pin the verified inode, not a later pathname open. Cryptography establishes
// trust once; kernel identity/change stamps detect accidental external writes
// afterwards (ctime as well as mtime, including nanoseconds). Each use checks
// both the descriptor and current directory entry before and after reading.
// This assumes an honest kernel/filesystem and no privileged timestamp forgery.
// It is not a cryptographic defence against a malicious kernel or root writer.
class TrustedArchive {
    int fd_ = -1;
    fs::path path_;
    struct stat identity_{};
    static bool equal(struct stat const & a, struct stat const & b) {
        return a.st_dev == b.st_dev && a.st_ino == b.st_ino && a.st_size == b.st_size &&
            a.st_mtim.tv_sec == b.st_mtim.tv_sec && a.st_mtim.tv_nsec == b.st_mtim.tv_nsec &&
            a.st_ctim.tv_sec == b.st_ctim.tv_sec && a.st_ctim.tv_nsec == b.st_ctim.tv_nsec &&
            S_ISREG(b.st_mode);
    }
public:
    explicit TrustedArchive(fs::path path) : path_{std::move(path)} {
        fd_ = open(path_.c_str(),O_RDONLY | O_NOFOLLOW | O_CLOEXEC | O_NONBLOCK);
        if (fd_ < 0) throw std::runtime_error("Cannot pin shard archive: " + path_.string());
        if (fstat(fd_,&identity_) || !S_ISREG(identity_.st_mode)) {
            close(fd_); fd_ = -1; throw std::runtime_error("Shard must be a regular file: " + path_.string());
        }
        try { check(); } catch (...) { close(fd_); fd_ = -1; throw; }
    }
    ~TrustedArchive() { if (fd_ >= 0) close(fd_); }
    TrustedArchive(TrustedArchive const &) = delete;
    TrustedArchive & operator=(TrustedArchive const &) = delete;
    fs::path pinned_path() const { return "/proc/self/fd/" + std::to_string(fd_); }
    std::size_t bytes() const { return static_cast<std::size_t>(identity_.st_size); }
    void check() const {
        struct stat descriptor{}, entry{};
        if (fstat(fd_,&descriptor) || lstat(path_.c_str(),&entry) ||
            !equal(identity_,descriptor) || !equal(identity_,entry))
            throw std::runtime_error("Shard changed while querying (identity/size/mtime/ctime): " + path_.string());
    }
};
