#include "config.hpp"
#include "logger.hpp"
#include "reference_index.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

namespace fs = std::filesystem;

namespace
{

struct NegativeSetConfig
{
    fs::path genome_file;
    fs::path output_dir{"."};
    fs::path output_fasta{"negative_set.fasta"};
    fs::path output_tsv{"negative_set.tsv"};
    fs::path log_file{"negative_set.log"};

    std::size_t kmer_size{15};
    std::size_t chromosome_count{1};
    std::size_t sequence_count{100};
    std::size_t min_length{80};
    std::size_t max_length{250};
    std::uint64_t seed{0};

    IBFConfig ibf{};
};

struct Chromosome
{
    std::string id;
    seqan3::dna5_vector sequence;
};

struct Sample
{
    std::size_t index{};
    std::size_t chromosome_index{};
    std::size_t start{};
    std::size_t end{};
    seqan3::dna5_vector original;
    seqan3::dna5_vector shuffled;
};

std::string trim(std::string_view text)
{
    auto const first = text.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos)
        return {};

    auto const last = text.find_last_not_of(" \t\r\n");
    return std::string{text.substr(first, last - first + 1)};
}

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

void validate_config(NegativeSetConfig const & cfg, fs::path const & config_path)
{
    if (cfg.genome_file.empty())
        throw std::runtime_error("Missing required config key 'genome_file' in " + config_path.string());

    if (cfg.kmer_size == 0 || cfg.kmer_size > 32)
        throw std::runtime_error("Config key 'kmer_size' must be in the range [1, 32].");

    if (cfg.chromosome_count == 0)
        throw std::runtime_error("Config key 'chromosome_count' must be at least 1.");

    if (cfg.sequence_count == 0)
        throw std::runtime_error("Config key 'sequence_count' must be at least 1.");

    if (cfg.min_length == 0)
        throw std::runtime_error("Config key 'min_length' must be at least 1.");

    if (cfg.max_length < cfg.min_length)
        throw std::runtime_error("Config key 'max_length' must be >= 'min_length'.");

    if (cfg.kmer_size > cfg.max_length)
        throw std::runtime_error("Config key 'kmer_size' must be <= 'max_length'.");

    if (cfg.ibf.hash_functions == 0 || cfg.ibf.hash_functions > 32)
        throw std::runtime_error("Config key 'ibf.hash_functions' must be in the range [1, 32].");

    if (cfg.ibf.fpr <= 0.0 || cfg.ibf.fpr > 0.5)
        throw std::runtime_error("Config key 'ibf.fpr' must be in the range (0, 0.5].");
}

NegativeSetConfig load_negative_set_config(fs::path const & config_path)
{
    std::ifstream in(config_path);
    if (!in)
        throw std::runtime_error("Failed to open config file: " + config_path.string());

    NegativeSetConfig cfg;
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
            if (current_section != "general" && current_section != "ibf")
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

        if (key == "genome_file" || key == "general.genome_file")
            cfg.genome_file = parse_string(key, value);
        else if (key == "output_dir" || key == "general.output_dir")
            cfg.output_dir = parse_string(key, value);
        else if (key == "output_fasta" || key == "general.output_fasta")
            cfg.output_fasta = parse_string(key, value);
        else if (key == "output_tsv" || key == "general.output_tsv")
            cfg.output_tsv = parse_string(key, value);
        else if (key == "log_file" || key == "general.log_file")
            cfg.log_file = parse_string(key, value);
        else if (key == "kmer_size" || key == "general.kmer_size")
            cfg.kmer_size = parse_size_t(key, value);
        else if (key == "chromosome_count" || key == "general.chromosome_count")
            cfg.chromosome_count = parse_size_t(key, value);
        else if (key == "sequence_count" || key == "general.sequence_count")
            cfg.sequence_count = parse_size_t(key, value);
        else if (key == "min_length" || key == "general.min_length")
            cfg.min_length = parse_size_t(key, value);
        else if (key == "max_length" || key == "general.max_length")
            cfg.max_length = parse_size_t(key, value);
        else if (key == "seed" || key == "general.seed")
            cfg.seed = parse_u64(key, value);
        else if (key == "ibf.hash_functions")
            cfg.ibf.hash_functions = parse_size_t(key, value);
        else if (key == "ibf.fpr")
            cfg.ibf.fpr = parse_double(key, value);
        else
            throw std::runtime_error("Unknown config key '" + key + "' in " + config_path.string());
    }

    validate_config(cfg, config_path);
    return cfg;
}

Config make_index_config(NegativeSetConfig const & cfg)
{
    Config index_cfg;
    index_cfg.index_method = IndexMethod::ibf;
    index_cfg.kmer_size = cfg.kmer_size;
    index_cfg.fragment_size = cfg.max_length;
    index_cfg.hit_threshold = 0;
    index_cfg.ibf = cfg.ibf;
    return index_cfg;
}

bool is_acgt(seqan3::dna5 const base)
{
    char const c = seqan3::to_char(base);
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

std::optional<seqan3::dna5_vector>
extract_clean_window(seqan3::dna5_vector const & seq, std::size_t start, std::size_t length)
{
    if (start > seq.size() || length > seq.size() - start)
        return std::nullopt;

    seqan3::dna5_vector window;
    window.reserve(length);

    for (std::size_t i = 0; i < length; ++i)
    {
        auto const base = seq[start + i];
        if (!is_acgt(base))
            return std::nullopt;
        window.push_back(base);
    }

    return window;
}

std::vector<Chromosome> read_chromosomes(fs::path const & genome_file, std::size_t chromosome_count)
{
    std::vector<Chromosome> chromosomes;
    chromosomes.reserve(chromosome_count);

    seqan3::sequence_file_input genome_in{genome_file};

    for (auto & record : genome_in)
    {
        Chromosome chromosome;
        chromosome.id = trim(record.id());
        chromosome.sequence = record.sequence();
        chromosomes.push_back(std::move(chromosome));

        if (chromosomes.size() == chromosome_count)
            break;
    }

    if (chromosomes.empty())
        throw std::runtime_error("No FASTA records found in genome file: " + genome_file.string());

    if (chromosomes.size() < chromosome_count)
        throw std::runtime_error("Requested " + std::to_string(chromosome_count) +
                                 " chromosome(s), but genome only contained " +
                                 std::to_string(chromosomes.size()) + " record(s).");

    return chromosomes;
}

std::vector<double> chromosome_weights(std::vector<Chromosome> const & chromosomes,
                                       std::size_t min_length)
{
    std::vector<double> weights;
    weights.reserve(chromosomes.size());

    for (auto const & chromosome : chromosomes)
        weights.push_back(chromosome.sequence.size() >= min_length
                              ? static_cast<double>(chromosome.sequence.size() - min_length + 1)
                              : 0.0);

    if (std::accumulate(weights.begin(), weights.end(), 0.0) == 0.0)
        throw std::runtime_error("None of the selected chromosomes is long enough for min_length.");

    return weights;
}

std::vector<Sample> sample_sequences(std::vector<Chromosome> const & chromosomes,
                                     NegativeSetConfig const & cfg,
                                     std::mt19937_64 & rng)
{
    auto weights = chromosome_weights(chromosomes, cfg.min_length);
    std::discrete_distribution<std::size_t> chromosome_dist(weights.begin(), weights.end()); //  produces random integers on the interval [0, n) https://en.cppreference.com/cpp/numeric/random/discrete_distribution
    std::uniform_int_distribution<std::size_t> length_dist(cfg.min_length, cfg.max_length);

    std::vector<Sample> samples;
    samples.reserve(cfg.sequence_count);

    std::size_t attempts = 0;
    std::size_t const max_attempts = std::max<std::size_t>(cfg.sequence_count * 1000, 10000);

    while (samples.size() < cfg.sequence_count && attempts < max_attempts)
    {
        ++attempts;

        std::size_t const chromosome_index = chromosome_dist(rng);
        auto const & chromosome = chromosomes[chromosome_index];
        if (chromosome.sequence.size() < cfg.min_length)
            continue;

        std::size_t const maximum_length = std::min(cfg.max_length, chromosome.sequence.size());
        std::size_t length = length_dist(rng);
        if (length > maximum_length)
            length = maximum_length;

        std::uniform_int_distribution<std::size_t> start_dist(0, chromosome.sequence.size() - length);
        std::size_t const start = start_dist(rng);
        auto window = extract_clean_window(chromosome.sequence, start, length);
        if (!window)
            continue;

        auto shuffled = *window;
        std::shuffle(shuffled.begin(), shuffled.end(), rng);

        samples.push_back(Sample{samples.size(),
                                 chromosome_index,
                                 start,
                                 start + length,
                                 *std::move(window),
                                 std::move(shuffled)});
    }

    if (samples.size() != cfg.sequence_count)
        throw std::runtime_error("Could only sample " + std::to_string(samples.size()) + " clean sequence(s) after " +
                                 std::to_string(attempts) + " attempt(s). Selected chromosomes may contain too many N bases.");

    return samples;
}

std::vector<seqan3::dna5_vector> originals_for_index(std::vector<Sample> const & samples)
{
    std::vector<seqan3::dna5_vector> originals;
    originals.reserve(samples.size());

    for (auto const & sample : samples)
        originals.push_back(sample.original);

    return originals;
}

std::uint64_t sum_counts(std::vector<std::size_t> const & counts)
{
    return std::accumulate(counts.begin(), counts.end(), std::uint64_t{0});
}

std::uint64_t count_nonzero(std::vector<std::size_t> const & counts)
{
    return static_cast<std::uint64_t>(std::count_if(counts.begin(), counts.end(),
                                                   [](std::size_t const count)
                                                   {
                                                       return count > 0;
                                                   }));
}

std::string output_id(Sample const & sample)
{
    return "negative_" + std::to_string(sample.index);
}

void write_outputs(std::vector<Chromosome> const & chromosomes,
                   std::vector<Sample> const & samples,
                   ReferenceIndexDna5 const & index,
                   NegativeSetConfig const & cfg)
{
    fs::create_directories(cfg.output_dir);

    auto fasta_path = cfg.output_dir / cfg.output_fasta;
    seqan3::sequence_file_output fasta_out{fasta_path};

    auto tsv_path = cfg.output_dir / cfg.output_tsv;
    std::ofstream tsv_out(tsv_path);
    if (!tsv_out)
        throw std::runtime_error("Failed to open TSV output file: " + tsv_path.string());

    tsv_out << "id\tchromosome\tstart\tend\tlength\ttotal_kmers\tmatched_kmer_positions\tmatch_count\tmatch_fraction\n";

    for (auto const & sample : samples)
    {
        auto counts = index.count_query_kmer_hits(sample.shuffled);
        std::uint64_t const total_kmers = sample.shuffled.size() >= cfg.kmer_size
                                              ? static_cast<std::uint64_t>(sample.shuffled.size() - cfg.kmer_size + 1)
                                              : 0;
        std::uint64_t const matched_positions = count_nonzero(counts);
        std::uint64_t const match_count = sum_counts(counts);
        double const match_fraction = total_kmers == 0 ? 0.0 : static_cast<double>(match_count) / total_kmers;

        auto const id = output_id(sample);
        fasta_out.emplace_back(sample.shuffled, id);

        tsv_out << id << '\t'
                << chromosomes[sample.chromosome_index].id << '\t'
                << sample.start << '\t'
                << sample.end << '\t'
                << sample.shuffled.size() << '\t'
                << total_kmers << '\t'
                << matched_positions << '\t'
                << match_count << '\t'
                << std::fixed << std::setprecision(6) << match_fraction << '\n';
    }

    Logger::info("Negative FASTA written to: " + fasta_path.string());
    Logger::info("Negative TSV written to: " + tsv_path.string());
}

} // namespace

int main(int argc, char ** argv)
{
    fs::path config_path = "negative_set/config.toml";

    try
    {
        if (argc > 2)
            throw std::runtime_error("Usage: ./negative_set [config.toml]");

        if (argc == 2)
            config_path = argv[1];

        auto cfg = load_negative_set_config(config_path);

        fs::create_directories(cfg.output_dir);
        Logger::init((cfg.output_dir / cfg.log_file).string());
        Logger::info("Starting negative_set.");
        Logger::info("Loaded config from: " + config_path.string());

        auto seed = cfg.seed;
        if (seed == 0)
            seed = static_cast<std::uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        std::mt19937_64 rng{seed};
        Logger::info("Using seed: " + std::to_string(seed));

        auto chromosomes = read_chromosomes(cfg.genome_file, cfg.chromosome_count);
        Logger::info("Loaded " + std::to_string(chromosomes.size()) + " chromosome record(s).");

        auto samples = sample_sequences(chromosomes, cfg, rng);
        Logger::info("Sampled and shuffled " + std::to_string(samples.size()) + " sequence(s).");

        auto index_cfg = make_index_config(cfg);
        auto originals = originals_for_index(samples);
        ReferenceIndexDna5 index{"negative_set_original_windows", originals, index_cfg};

        write_outputs(chromosomes, samples, index, cfg);
        Logger::info("negative_set finished successfully.");
    }
    catch (std::exception const & e)
    {
        Logger::error(std::string{"Fatal error: "} + e.what());
        std::cerr << "Fatal error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
