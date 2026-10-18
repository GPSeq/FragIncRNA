#pragma once

#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <atomic>
#include <csignal>
#include <algorithm>
#include <filesystem>
#include <map>
#include <vector>
#include <unistd.h>
#include <sys/resource.h>
#if defined(__linux__)
#include <sched.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace fraglncrna {
inline std::atomic<int> requested_signal{0};
static_assert(std::atomic<int>::is_always_lock_free);
inline void request_stop(int signal) { requested_signal.store(signal,std::memory_order_relaxed); }
inline void check_stop() {
    auto signal = requested_signal.load(std::memory_order_relaxed);
    if (signal) throw std::runtime_error("Run interrupted by signal " + std::to_string(signal));
}
inline std::string log_escape(std::string const & text) {
    std::ostringstream out;
    for (unsigned char c : text) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << unsigned(c);
        else out << c;
    }
    return out.str();
}

// Coordinator/worker safe, bounded, flushed line-at-a-time logs. Scientific
// output streams never pass through this class. No background writer to drain.
class ProductionLog {
    std::mutex mutex_;
    std::mutex progress_mutex_;
    std::ofstream file_, json_;
    std::chrono::steady_clock::time_point start_ = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point last_progress_ = start_;
    int threshold_ = 2;
    bool console_ = true;
    double interval_ = 10;
    std::string run_id_;
    std::map<std::string, std::chrono::steady_clock::time_point> stages_;
public:
    std::size_t warnings = 0, errors = 0;
    std::string const & run_id() const { return run_id_; }
    void open(std::filesystem::path const & path, std::string const & level, bool console, double interval) {
        console_ = console; interval_ = interval;
        std::vector<std::string> levels{"trace", "debug", "info", "warning", "error", "off"};
        auto it = std::find(levels.begin(), levels.end(), level);
        if (it == levels.end()) throw std::runtime_error("Invalid --log-level (trace/debug/info/warning/error/off)");
        threshold_ = int(it - levels.begin());
        run_id_ = std::to_string(getpid()) + "-" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count());
        if (threshold_ == 5) return;
        std::filesystem::create_directories(path.parent_path().empty() ? "." : path.parent_path());
        file_.open(path); json_.open(path.string() + ".jsonl");
        if (!file_ || !json_) throw std::runtime_error("Cannot write log " + path.string());
    }
    void event(std::string const & phase, std::string const & event, std::string const & message, int severity = 2) {
        std::lock_guard lock(mutex_);
        if (severity == 3) ++warnings;
        if (severity == 4) ++errors;
        auto steady = std::chrono::steady_clock::now();
        if (event == "stage_started") stages_[phase] = steady;
        if (severity < threshold_ || threshold_ == 5) return;
        auto now = std::chrono::system_clock::now();
        auto stamp = std::chrono::system_clock::to_time_t(now);
        std::tm tm{}; gmtime_r(&stamp, &tm);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count() % 1000;
        std::ostringstream time;
        time << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S") << '.' << std::setw(3) << std::setfill('0') << ms << 'Z';
        double total = std::chrono::duration<double>(steady - start_).count();
        double phase_seconds = stages_.contains(phase) ? std::chrono::duration<double>(steady - stages_[phase]).count() : 0;
        static char const * labels[]{"TRACE", "DEBUG", "INFO", "WARN", "ERROR"};
        std::ostringstream line;
        line << time.str() << ' ' << labels[severity] << " [run=" << run_id_ << "] [" << phase << "] "
             << event << " elapsed=" << total << " phase_seconds=" << phase_seconds << ' ' << message << '\n';
        if (console_) std::cerr << line.str();
        if (file_.is_open()) {
            file_ << line.str(); file_.flush();
            json_ << "{\"timestamp\":\"" << time.str() << "\",\"run_id\":\"" << run_id_
                  << "\",\"level\":\"" << labels[severity] << "\",\"phase\":\"" << phase
                  << "\",\"event\":\"" << event << "\",\"elapsed_seconds\":" << total
                  << ",\"phase_seconds\":" << phase_seconds << ",\"message\":\"" << log_escape(message) << "\"}\n";
            json_.flush();
            if (!file_ || !json_) throw std::runtime_error("Failed writing production logs; check free disk space");
        }
    }
    void progress(std::string const & phase, std::string const & message, bool force = false) {
        {
            std::lock_guard lock(progress_mutex_);
            auto now = std::chrono::steady_clock::now();
            if (!force && std::chrono::duration<double>(now - last_progress_).count() < interval_) return;
            last_progress_ = now;
        }
        event(phase, "stage_progress", message);
    }
};

inline ProductionLog production_log;

inline int available_workers() {
    int cpus = std::max(1u, std::thread::hardware_concurrency());
#if defined(__linux__)
    cpu_set_t mask; CPU_ZERO(&mask);
    if (sched_getaffinity(0, sizeof(mask), &mask) == 0) cpus = std::max(1, CPU_COUNT(&mask));
    std::ifstream quota("/sys/fs/cgroup/cpu.max");
    std::string amount; long period = 0;
    if (quota >> amount >> period && amount != "max" && period > 0)
        cpus = std::min(cpus, std::max(1, int((std::stoll(amount) + period - 1) / period)));
#endif
    if (auto slurm = std::getenv("SLURM_CPUS_PER_TASK")) cpus = std::min(cpus, std::max(1, std::atoi(slurm)));
    return cpus;
}

template <typename Function>
void parallel_records(std::size_t count, [[maybe_unused]] int workers, Function function) {
    std::exception_ptr error;
    std::mutex error_mutex;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(workers) if(workers > 1)
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(count); ++i) {
        try { check_stop(); function(static_cast<std::size_t>(i)); }
        catch (...) { std::lock_guard lock(error_mutex); if (!error) error = std::current_exception(); }
    }
    if (error) std::rethrow_exception(error);
}
} // namespace fraglncrna
