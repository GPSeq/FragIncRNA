#include "logger.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

std::ofstream Logger::file_;
std::mutex Logger::mutex_;

namespace
{

/*
* @fn current_timestamp
* @brief Formats the current local time for log entries.
* @signature std::string current_timestamp();
* @param None.
* @throws None.
* @return Timestamp string formatted as YYYY-MM-DD HH:MM:SS.
*/
std::string current_timestamp()
{
    using namespace std::chrono;
    auto now = system_clock::now();
    std::time_t t = system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

/*
* @fn log_impl
* @brief Writes a timestamped log line to the configured log file.
* @signature void log_impl(std::string const & level, std::string const & msg, std::ofstream & file, std::mutex & mtx);
* @param level: log severity label.
* @param msg: message text to log.
* @param file: output log file stream.
* @param mtx: mutex protecting the log stream.
* @throws None.
* @return None.
*/
void log_impl(std::string const & level, std::string const & msg,
              std::ofstream & file, std::mutex & mtx)
{
    std::lock_guard<std::mutex> lock{mtx};
    std::string line = "[" + level + "] " + current_timestamp() + " " + msg;
    // no stdout here: keep stdout clean for progress bar only
    if (file.is_open())
        file << line << '\n';
}


} // namespace

/*
* @fn init
* @brief Opens the log file and resets any previously open log stream.
* @signature void Logger::init(std::string const & path);
* @param path: path to the log file to create or truncate.
* @throws None.
* @return None.
*/
void Logger::init(std::string const & path)
{
    std::lock_guard<std::mutex> lock{mutex_};
    if (file_.is_open())
        file_.close();
    file_.open(path, std::ios::out | std::ios::trunc);
}

/*
* @fn info
* @brief Writes an informational message to the log file.
* @signature void Logger::info(std::string const & msg);
* @param msg: message text to log.
* @throws None.
* @return None.
*/
void Logger::info(std::string const & msg)
{
    log_impl("INFO", msg, file_, mutex_);
}

/*
* @fn warn
* @brief Writes a warning message to the log file.
* @signature void Logger::warn(std::string const & msg);
* @param msg: message text to log.
* @throws None.
* @return None.
*/
void Logger::warn(std::string const & msg)
{
    log_impl("WARN", msg, file_, mutex_);
}

/*
* @fn error
* @brief Writes an error message to the log file.
* @signature void Logger::error(std::string const & msg);
* @param msg: message text to log.
* @throws None.
* @return None.
*/
void Logger::error(std::string const & msg)
{
    log_impl("ERROR", msg, file_, mutex_);
}

/*
* @fn print_stdout
* @brief Writes a progress or status message to standard output.
* @signature void Logger::print_stdout(std::string const & msg, bool newline);
* @param msg: message text to print.
* @param newline: whether to append a newline after the message.
* @throws None.
* @return None.
*/
void Logger::print_stdout(std::string const & msg, bool newline)
{
    std::lock_guard<std::mutex> lock{mutex_};
    std::cout << msg;
    if (newline)
        std::cout << '\n';
    std::cout << std::flush;
}
