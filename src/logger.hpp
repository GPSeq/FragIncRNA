#pragma once

/**
 * @file logger.hpp
 * @brief Thread-safe text logger of the legacy executable `lncrna_mers`.
 */
#include <fstream>
#include <mutex>
#include <string>

class Logger
{
public:
    /**
     * @brief Opens the log file and resets any previously open log stream.
     * @param path path to the log file to create or truncate.
     */
    static void init(std::string const & path);

    /**
     * @brief Writes an informational message to the log file.
     * @param msg message text to log.
     */
    static void info(std::string const & msg);

    /**
     * @brief Writes a warning message to the log file.
     * @param msg message text to log.
     */
    static void warn(std::string const & msg);

    /**
     * @brief Writes an error message to the log file.
     * @param msg message text to log.
     */
    static void error(std::string const & msg);

    /**
     * @brief Writes a progress or status message to standard output.
     * @param msg message text to print.
     * @param newline whether to append a newline after the message.
     */
    static void print_stdout(std::string const & msg, bool newline = false);

private:
    static std::ofstream file_;
    static std::mutex mutex_;
};
