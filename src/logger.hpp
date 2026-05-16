#pragma once

#include <fstream>
#include <mutex>
#include <string>

class Logger
{
public:
    /*
    * @fn init
    * @brief Opens the log file and resets any previously open log stream.
    * @signature static void init(std::string const & path);
    * @param path: path to the log file to create or truncate.
    * @throws None.
    * @return None.
    */
    static void init(std::string const & path);

    /*
    * @fn info
    * @brief Writes an informational message to the log file.
    * @signature static void info(std::string const & msg);
    * @param msg: message text to log.
    * @throws None.
    * @return None.
    */
    static void info(std::string const & msg);

    /*
    * @fn warn
    * @brief Writes a warning message to the log file.
    * @signature static void warn(std::string const & msg);
    * @param msg: message text to log.
    * @throws None.
    * @return None.
    */
    static void warn(std::string const & msg);

    /*
    * @fn error
    * @brief Writes an error message to the log file.
    * @signature static void error(std::string const & msg);
    * @param msg: message text to log.
    * @throws None.
    * @return None.
    */
    static void error(std::string const & msg);

    /*
    * @fn print_stdout
    * @brief Writes a progress or status message to standard output.
    * @signature static void print_stdout(std::string const & msg, bool newline = false);
    * @param msg: message text to print.
    * @param newline: whether to append a newline after the message.
    * @throws None.
    * @return None.
    */
    static void print_stdout(std::string const & msg, bool newline = false);

private:
    static std::ofstream file_;
    static std::mutex mutex_;
};
