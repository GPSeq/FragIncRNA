#include "config.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace
{

/*
* @fn trim
* @brief Removes leading and trailing whitespace from a string view.
* @signature std::string trim(std::string_view text);
* @param text: text to trim.
* @throws None.
* @return Trimmed string.
*/
std::string trim(std::string_view text)
{
    auto const first = text.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos)
        return {};

    auto const last = text.find_last_not_of(" \t\r\n");
    return std::string{text.substr(first, last - first + 1)};
}

/*
* @fn strip_comment
* @brief Removes TOML comments from a line while preserving hash characters inside quoted strings.
* @signature std::string strip_comment(std::string_view line);
* @param line: TOML line to process.
* @throws None.
* @return Line content before an unquoted comment marker.
*/
std::string strip_comment(std::string_view line)
{
    bool in_quotes = false;

    for (std::size_t i = 0; i < line.size(); ++i)
    {
        if (line[i] == '"' && (i == 0 || line[i - 1] != '\\'))
            in_quotes = !in_quotes;

        if (!in_quotes && line[i] == '#')
            return std::string{line.substr(0, i)};
    }

    return std::string{line};
}

/*
* @fn parse_string
* @brief Parses a TOML double-quoted string value and decodes supported escape sequences.
* @signature std::string parse_string(std::string const & key, std::string const & value);
* @param key: configuration key used in error messages.
* @param value: raw TOML value to parse.
* @throws std::runtime_error when the value is not a valid supported TOML string.
* @return Parsed string value.
*/
std::string parse_string(std::string const & key, std::string const & value)
{
    if (value.size() < 2 || value.front() != '"' || value.back() != '"')
        throw std::runtime_error("Config key '" + key + "' must be a TOML string enclosed in double quotes.");

    std::string out;
    out.reserve(value.size() - 2);

    for (std::size_t i = 1; i + 1 < value.size(); ++i)
    {
        if (value[i] == '\\' && i + 1 < value.size() - 1)
        {
            ++i;
            switch (value[i])
            {
                case '\\': out.push_back('\\'); break;
                case '"': out.push_back('"'); break;
                case 'n': out.push_back('\n'); break;
                case 't': out.push_back('\t'); break;
                default:
                    throw std::runtime_error("Unsupported escape sequence in config key '" + key + "'.");
            }
            continue;
        }

        out.push_back(value[i]);
    }

    return out;
}

/*
* @fn parse_bool
* @brief Parses a TOML boolean value.
* @signature bool parse_bool(std::string const & key, std::string const & value);
* @param key: configuration key used in error messages.
* @param value: raw TOML value to parse.
* @throws std::runtime_error when the value is not true or false.
* @return Parsed boolean value.
*/
bool parse_bool(std::string const & key, std::string const & value)
{
    if (value == "true")
        return true;
    if (value == "false")
        return false;

    throw std::runtime_error("Config key '" + key + "' must be true or false.");
}

/*
* @fn parse_size_t
* @brief Parses a non-negative integer into std::size_t.
* @signature std::size_t parse_size_t(std::string const & key, std::string const & value);
* @param key: configuration key used in error messages.
* @param value: raw TOML value to parse.
* @throws std::runtime_error when the value is not a valid non-negative integer.
* @return Parsed std::size_t value.
*/
std::size_t parse_size_t(std::string const & key, std::string const & value)
{
    try
    {
        std::size_t pos{};
        auto const parsed = std::stoull(value, &pos);
        if (pos != value.size())
            throw std::runtime_error("");
        return static_cast<std::size_t>(parsed);
    }
    catch (...)
    {
        throw std::runtime_error("Config key '" + key + "' must be a non-negative integer.");
    }
}

/*
* @fn parse_u64
* @brief Parses a non-negative integer into std::uint64_t.
* @signature std::uint64_t parse_u64(std::string const & key, std::string const & value);
* @param key: configuration key used in error messages.
* @param value: raw TOML value to parse.
* @throws std::runtime_error when the value is not a valid non-negative integer.
* @return Parsed std::uint64_t value.
*/
std::uint64_t parse_u64(std::string const & key, std::string const & value)
{
    try
    {
        std::size_t pos{};
        auto const parsed = std::stoull(value, &pos);
        if (pos != value.size())
            throw std::runtime_error("");
        return parsed;
    }
    catch (...)
    {
        throw std::runtime_error("Config key '" + key + "' must be a non-negative integer.");
    }
}

/*
* @fn parse_double
* @brief Parses a floating-point TOML value.
* @signature double parse_double(std::string const & key, std::string const & value);
* @param key: configuration key used in error messages.
* @param value: raw TOML value to parse.
* @throws std::runtime_error when the value is not a valid floating-point number.
* @return Parsed double value.
*/
double parse_double(std::string const & key, std::string const & value)
{
    try
    {
        std::size_t pos{};
        auto const parsed = std::stod(value, &pos);
        if (pos != value.size())
            throw std::runtime_error("");
        return parsed;
    }
    catch (...)
    {
        throw std::runtime_error("Config key '" + key + "' must be a floating-point number.");
    }
}

/*
* @fn parse_index_method
* @brief Parses the configured index method.
* @signature IndexMethod parse_index_method(std::string const & value);
* @param value: parsed TOML string containing the index method name.
* @throws std::runtime_error when the value is not a supported index method.
* @return Parsed index method.
*/
IndexMethod parse_index_method(std::string const & value)
{
    if (value == "ibf")
        return IndexMethod::ibf;
    if (value == "hibf")
        return IndexMethod::hibf;

    throw std::runtime_error("Config key 'general.index_method' must be \"ibf\" or \"hibf\".");
}

/*
* @fn validate_config
* @brief Validates required configuration values and option ranges.
* @signature void validate_config(Config const & cfg, std::filesystem::path const & config_path);
* @param cfg: configuration to validate.
* @param config_path: source configuration path used in error messages.
* @throws std::runtime_error when any required value is missing or out of range.
* @return None.
*/
void validate_config(Config const & cfg, std::filesystem::path const & config_path)
{
    if (cfg.ref_dir.empty())
        throw std::runtime_error("Missing required config key 'ref_dir' in " + config_path.string());

    if (cfg.query_file.empty())
        throw std::runtime_error("Missing required config key 'query_file' in " + config_path.string());

    if (cfg.fragment_size < 4)
        throw std::runtime_error("Config key 'fragment_size' must be at least 4.");

    if (cfg.kmer_size == 0 || cfg.kmer_size > 32)
        throw std::runtime_error("Config key 'kmer_size' must be in the range [1, 32].");

    if (cfg.index_method == IndexMethod::ibf)
    {
        if (cfg.ibf.hash_functions == 0 || cfg.ibf.hash_functions > 32)
            throw std::runtime_error("Config key 'ibf.hash_functions' must be in the range [1, 32].");

        if (cfg.ibf.fpr <= 0.0 || cfg.ibf.fpr > 0.5)
            throw std::runtime_error("Config key 'ibf.fpr' must be in the range (0, 0.5].");
    }
    else
    {
        if (cfg.hibf.hash_functions == 0 || cfg.hibf.hash_functions > 5)
            throw std::runtime_error("Config key 'hibf.hash_functions' must be in the range [1, 5].");

        if (cfg.hibf.maximum_fpr <= 0.0 || cfg.hibf.maximum_fpr >= 1.0)
            throw std::runtime_error("Config key 'hibf.maximum_fpr' must be in the range (0, 1).");

        if (cfg.hibf.relaxed_fpr <= 0.0 || cfg.hibf.relaxed_fpr >= 1.0)
            throw std::runtime_error("Config key 'hibf.relaxed_fpr' must be in the range (0, 1).");

        if (cfg.hibf.relaxed_fpr < cfg.hibf.maximum_fpr)
            throw std::runtime_error("Config key 'hibf.relaxed_fpr' must be >= 'hibf.maximum_fpr'.");

        if (cfg.hibf.sketch_bits < 5 || cfg.hibf.sketch_bits > 32)
            throw std::runtime_error("Config key 'hibf.sketch_bits' must be in the range [5, 32].");

        if (cfg.hibf.empty_bin_fraction < 0.0 || cfg.hibf.empty_bin_fraction >= 1.0)
            throw std::runtime_error("Config key 'hibf.empty_bin_fraction' must be in the range [0.0, 1.0).");

        if (cfg.hibf.alpha <= 0.0)
            throw std::runtime_error("Config key 'hibf.alpha' must be positive.");

        if (cfg.hibf.max_rearrangement_ratio < 0.0 || cfg.hibf.max_rearrangement_ratio > 1.0)
            throw std::runtime_error("Config key 'hibf.max_rearrangement_ratio' must be in the range [0.0, 1.0].");
    }
}

} // namespace

/*
* @fn load_config_from_toml
* @brief Loads and validates application configuration from a TOML file.
* @signature Config load_config_from_toml(std::filesystem::path const & config_path);
* @param config_path: path to the TOML configuration file.
* @throws std::runtime_error when the file cannot be opened or contains invalid configuration.
* @return Parsed and validated configuration.
*/
Config load_config_from_toml(std::filesystem::path const & config_path)
{
    std::ifstream in(config_path);
    if (!in)
        throw std::runtime_error("Failed to open config file: " + config_path.string());

    Config cfg;
    std::unordered_set<std::string> seen_keys;
    std::string current_section;

    std::string raw_line;
    std::size_t line_number = 0;
    while (std::getline(in, raw_line))
    {
        ++line_number;
        auto line = trim(strip_comment(raw_line));
        if (line.empty())
            continue;

        if (line.front() == '[')
        {
            if (line.back() != ']')
                throw std::runtime_error("Malformed TOML section header at " + config_path.string() +
                                         ":" + std::to_string(line_number));

            current_section = trim(std::string_view{line}.substr(1, line.size() - 2));
            if (current_section != "general" && current_section != "ibf" && current_section != "hibf")
                throw std::runtime_error("Unsupported TOML section [" + current_section + "] in " +
                                         config_path.string() + ":" + std::to_string(line_number));
            continue;
        }

        auto const eq_pos = line.find('=');
        if (eq_pos == std::string::npos)
            throw std::runtime_error("Expected key = value at " + config_path.string() +
                                     ":" + std::to_string(line_number));

        auto raw_key = trim(std::string_view{line}.substr(0, eq_pos));
        auto value = trim(std::string_view{line}.substr(eq_pos + 1));

        if (raw_key.empty() || value.empty())
            throw std::runtime_error("Invalid key/value pair at " + config_path.string() +
                                     ":" + std::to_string(line_number));

        auto key = current_section.empty() ? raw_key : current_section + "." + raw_key;

        if (!seen_keys.insert(key).second)
            throw std::runtime_error("Duplicate config key '" + key + "' in " + config_path.string());

        if (key == "ref_dir" || key == "general.ref_dir")
            cfg.ref_dir = parse_string(key, value);
        else if (key == "query_file" || key == "general.query_file")
            cfg.query_file = parse_string(key, value);
        else if (key == "output_dir" || key == "general.output_dir")
            cfg.output_dir = parse_string(key, value);
        else if (key == "output_file" || key == "general.output_file")
            cfg.output_file = parse_string(key, value);
        else if (key == "log_file" || key == "general.log_file")
            cfg.log_file = parse_string(key, value);
        else if (key == "index_method" || key == "general.index_method")
            cfg.index_method = parse_index_method(parse_string(key, value));
        else if (key == "fragment_size" || key == "general.fragment_size")
            cfg.fragment_size = parse_size_t(key, value);
        else if (key == "kmer_size" || key == "general.kmer_size")
            cfg.kmer_size = parse_size_t(key, value);
        else if (key == "hit_threshold" || key == "general.hit_threshold")
            cfg.hit_threshold = parse_u64(key, value);
        else if (key == "store_fragments" || key == "general.store_fragments")
            cfg.store_fragments = parse_bool(key, value);
        else if (key == "store_index" || key == "general.store_index")
            cfg.store_index = parse_bool(key, value);
        else if (key == "cleanup_index" || key == "general.cleanup_index")
            cfg.cleanup_index = parse_bool(key, value);
        else if (key == "single_results_writer" || key == "general.single_results_writer")
            cfg.single_results_writer = parse_bool(key, value);
        else if (key == "threads" || key == "general.threads")
            cfg.threads = parse_size_t(key, value);
        else if (key == "ibf.hash_functions")
            cfg.ibf.hash_functions = parse_size_t(key, value);
        else if (key == "ibf.fpr")
            cfg.ibf.fpr = parse_double(key, value);
        else if (key == "hibf.hash_functions")
            cfg.hibf.hash_functions = parse_size_t(key, value);
        else if (key == "hibf.maximum_fpr")
            cfg.hibf.maximum_fpr = parse_double(key, value);
        else if (key == "hibf.relaxed_fpr")
            cfg.hibf.relaxed_fpr = parse_double(key, value);
        else if (key == "hibf.threads")
            cfg.hibf.threads = parse_size_t(key, value);
        else if (key == "hibf.sketch_bits")
            cfg.hibf.sketch_bits = parse_size_t(key, value);
        else if (key == "hibf.tmax")
            cfg.hibf.tmax = parse_size_t(key, value);
        else if (key == "hibf.empty_bin_fraction")
            cfg.hibf.empty_bin_fraction = parse_double(key, value);
        else if (key == "hibf.alpha")
            cfg.hibf.alpha = parse_double(key, value);
        else if (key == "hibf.max_rearrangement_ratio")
            cfg.hibf.max_rearrangement_ratio = parse_double(key, value);
        else if (key == "hibf.disable_estimate_union")
            cfg.hibf.disable_estimate_union = parse_bool(key, value);
        else if (key == "hibf.disable_rearrangement")
            cfg.hibf.disable_rearrangement = parse_bool(key, value);
        else
            throw std::runtime_error("Unknown config key '" + key + "' in " + config_path.string());
    }

    validate_config(cfg, config_path);
    return cfg;
}
