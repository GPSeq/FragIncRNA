// Version 2 analysis path. The historical executable and its TSV remain available.
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include "sha256.hpp"
#include "kmer_windows.hpp"
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <zlib.h>

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <map>
#include <limits>
#include <numeric>
#include <optional>
#include <ranges>
#include <tuple>
#include <stdexcept>
#include <sstream>
#include <string>
#include <string_view>
#include <sys/resource.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <unistd.h>

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

namespace
{
struct Options
{
    fs::path reference, queries, output, metrics, config_path;
    fs::path positions, calibration_bins, fragments;
    fs::path reference_repeat_bed, query_repeat_bed, repeat_classes, reference_hit_classes, exon_blocks;
    fs::path loci_output, expected_loci_bed, expected_loci_output;
    fs::path save_index, load_index;
    fs::path threshold_table;
    std::size_t k = 15, fragment_size = 8000, hash_functions = 3, batch_windows = 100000;
    std::size_t local_window = 0, max_position_rows = 1000000, max_calibration_pairs = 10000000;
    std::size_t max_ibf_bin_pairs = 10000000;
    std::size_t max_repeat_intervals = 5000000, max_exon_blocks = 1000000, max_query_bases = 50000000;
    std::size_t locus_gap = 8000, max_locus_span = 1000000, max_locus_rows = 100000, max_expected_loci = 1000000;
    bool local_window_given = false, genome_fpr_given = false;
    double fpr = 0.01, basic = 0.75, strict = 0.85;
    double genome_fpr = 0.01, requested_per_bin_fpr = 0.01;
    std::string strand = "canonical", ambiguity = "skip", repeat = "include", mode = "validate";
    std::string mask_target = "both";
    std::string evaluation_control_set, training_control_set, threshold_source = "cli";
};

std::string reference_id_from_path(fs::path const & path)
{
    std::string id = path.stem().string();
    return path.extension() == ".gz" ? fs::path{id}.stem().string() : id;
}

std::vector<std::string> split_tabs(std::string const & line)
{
    std::vector<std::string> fields;
    std::size_t begin = 0, end;
    while ((end = line.find('\t', begin)) != std::string::npos)
    {
        fields.push_back(line.substr(begin, end - begin));
        begin = end + 1;
    }
    fields.push_back(line.substr(begin));
    return fields;
}

void load_thresholds(Options & o)
{
    if (o.threshold_table.empty()) return;
    if (o.evaluation_control_set.empty())
        throw std::runtime_error("--evaluation-control-set is required with --threshold-table");
    std::ifstream input(o.threshold_table);
    if (!input) throw std::runtime_error("Cannot open threshold table " + o.threshold_table.string());
    std::string line;
    if (!std::getline(input, line))
        throw std::runtime_error("Empty threshold table");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line != "reference_id\tk\tbasic_threshold\tstrict_threshold\ttraining_control_set")
        throw std::runtime_error("Threshold table header must be reference_id, k, basic_threshold, "
                                 "strict_threshold, training_control_set (TSV)");
    std::string ref_id = reference_id_from_path(o.reference);
    bool selected = false;
    while (std::getline(input, line))
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        auto fields = split_tabs(line);
        if (fields.size() != 5) throw std::runtime_error("Invalid threshold table row");
        if (fields[0] != ref_id || std::stoull(fields[1]) != o.k) continue;
        if (selected) throw std::runtime_error("Duplicate threshold table row for reference/k");
        o.basic = std::stod(fields[2]);
        o.strict = std::stod(fields[3]);
        o.training_control_set = fields[4];
        if (o.training_control_set.empty() || o.training_control_set == o.evaluation_control_set)
            throw std::runtime_error("Training and evaluation control sets must be nonempty and different");
        selected = true;
    }
    if (!selected) throw std::runtime_error("No threshold table row for reference/k");
    if (!(0 <= o.basic && o.basic <= o.strict && o.strict <= 1))
        throw std::runtime_error("Invalid threshold values in threshold table");
    o.threshold_source = "table";
}

std::string json_string(std::string const & value)
{
    std::string out = "\"";
    for (char c : value)
    {
        if (c == '"' || c == '\\') out += '\\';
        if (c == '\n') out += "\\n";
        else if (c == '\t') out += "\\t";
        else if (static_cast<unsigned char>(c) < 32) out += '?';
        else out += c;
    }
    return out + '"';
}

uint64_t fnv_update(uint64_t hash, unsigned char c)
{
    return (hash ^ c) * 1099511628211ULL;
}

std::string hex64(uint64_t value)
{
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16) << value;
    return out.str();
}

std::string sha256_file(fs::path const & path, std::uint64_t offset = 0)
{
    if (path.empty()) return {};
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot hash " + path.string());
    if (offset)
    {
        if (offset > fs::file_size(path)) throw std::runtime_error("Invalid hash offset for " + path.string());
        input.seekg(static_cast<std::streamoff>(offset));
    }
    fraglncrna::Sha256 sha;
    std::array<char, 65536> buffer{};
    while (input.read(buffer.data(), buffer.size()) || input.gcount())
        for (std::streamsize i = 0; i < input.gcount(); ++i)
            sha.update(static_cast<unsigned char>(buffer[static_cast<std::size_t>(i)]));
    if (!input.eof()) throw std::runtime_error("Failed reading " + path.string());
    return sha.hex_digest();
}

std::string sha256_text(std::string_view value)
{
    fraglncrna::Sha256 sha;
    sha.update(value);
    return sha.hex_digest();
}

void sha256_u64(fraglncrna::Sha256 & sha, uint64_t value)
{
    for (unsigned i = 0; i < 8; ++i) sha.update(static_cast<unsigned char>(value >> (i * 8)));
}

void sha256_record_start(fraglncrna::Sha256 & sha, std::string const & id)
{
    sha.update(0xff);
    sha256_u64(sha, id.size());
    sha.update(id);
}

void sha256_record_end(fraglncrna::Sha256 & sha, uint64_t sequence_length)
{
    sha.update(0xfe);
    sha256_u64(sha, sequence_length);
}

std::string trim_config(std::string_view text)
{
    auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos) return {};
    auto last = text.find_last_not_of(" \t\r\n");
    return std::string{text.substr(first, last - first + 1)};
}

std::string strip_config_comment(std::string_view line)
{
    char quote = 0;
    bool escaped = false;
    for (std::size_t i = 0; i < line.size(); ++i)
    {
        char c = line[i];
        if (quote == '"' && escaped) { escaped = false; continue; }
        if (quote == '"' && c == '\\') { escaped = true; continue; }
        if (c == quote) quote = 0;
        else if (!quote && (c == '"' || c == '\'')) quote = c;
        else if (!quote && c == '#') return std::string{line.substr(0, i)};
    }
    return std::string{line};
}

std::string parse_config_string(std::string const & key, std::string const & raw)
{
    if (raw.size() < 2 || (raw.front() != '"' && raw.front() != '\'') || raw.back() != raw.front())
        throw std::runtime_error("TOML key '" + key + "' requires a quoted string");
    if (raw.front() == '\'') return raw.substr(1, raw.size() - 2);
    std::string result;
    for (std::size_t i = 1; i + 1 < raw.size(); ++i)
    {
        char c = raw[i];
        if (c == '"') throw std::runtime_error("Unescaped quote in TOML key '" + key + "'");
        if (c != '\\') { result += c; continue; }
        if (++i + 1 >= raw.size()) throw std::runtime_error("Invalid TOML escape in '" + key + "'");
        switch (raw[i])
        {
            case '\\': result += '\\'; break;
            case '"': result += '"'; break;
            case 'n': result += '\n'; break;
            case 'r': result += '\r'; break;
            case 't': result += '\t'; break;
            default: throw std::runtime_error("Unsupported TOML escape in '" + key + "'");
        }
    }
    return result;
}

std::vector<std::string> load_v2_config_args(fs::path const & path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open TOML config " + path.string());
    std::vector<std::string> args;
    std::unordered_set<std::string> seen;
    std::string section, raw;
    std::size_t line_number = 0;
    auto is_integer = [](std::string const & key)
    {
        return key == "k" || key == "fragment_size" || key == "hash_functions" ||
               key == "batch_windows" || key == "local_window" || key == "max_position_rows" ||
               key == "max_calibration_pairs" || key == "max_ibf_bin_pairs" ||
               key == "max_repeat_intervals" || key == "max_exon_blocks" ||
               key == "max_query_bases" || key == "locus_gap" || key == "max_locus_span" ||
               key == "max_locus_rows" || key == "max_expected_loci";
    };
    auto is_real = [](std::string const & key)
    {
        return key == "fpr" || key == "genome_fpr" || key == "basic_threshold" ||
               key == "strict_threshold";
    };
    auto is_string = [](std::string const & key)
    {
        return key == "reference" || key == "queries" || key == "output" || key == "metrics" ||
               key == "positions" || key == "calibration_bins" || key == "fragments" ||
               key == "reference_repeat_bed" || key == "query_repeat_bed" ||
               key == "repeat_classes" || key == "reference_hit_classes" || key == "exon_blocks" ||
               key == "discover_loci" || key == "expected_loci_bed" ||
               key == "expected_loci_output" || key == "save_index" || key == "load_index" ||
               key == "mask_target" || key == "threshold_table" ||
               key == "evaluation_control_set" || key == "strand" || key == "ambiguity" ||
               key == "repeat" || key == "mode";
    };
    while (std::getline(input, raw))
    {
        ++line_number;
        std::string line = trim_config(strip_config_comment(raw));
        if (line.empty()) continue;
        if (line.front() == '[')
        {
            if (line != "[v2]" || !section.empty())
                throw std::runtime_error("Expected one [v2] TOML section at " + path.string() + ":" +
                                         std::to_string(line_number));
            section = "v2";
            continue;
        }
        if (section != "v2")
            throw std::runtime_error("TOML keys must be inside [v2] at " + path.string() + ":" +
                                     std::to_string(line_number));
        auto equal = line.find('=');
        if (equal == std::string::npos)
            throw std::runtime_error("Expected key = value at " + path.string() + ":" +
                                     std::to_string(line_number));
        std::string key = trim_config(std::string_view{line}.substr(0, equal));
        std::string value = trim_config(std::string_view{line}.substr(equal + 1));
        if (!seen.insert(key).second)
            throw std::runtime_error("Duplicate TOML key '" + key + "'");
        if (!is_string(key) && !is_integer(key) && !is_real(key))
            throw std::runtime_error("Unknown v2 TOML key '" + key + "'");
        if (is_string(key)) value = parse_config_string(key, value);
        else
        {
            if (value.empty() || value.front() == '-' || value.front() == '"' || value.front() == '\'')
                throw std::runtime_error("TOML key '" + key + "' requires a nonnegative number");
            std::string number;
            for (std::size_t j = 0; j < value.size(); ++j)
            {
                char c = value[j];
                if (c == '_')
                {
                    if (j == 0 || j + 1 == value.size() ||
                        !std::isdigit(static_cast<unsigned char>(value[j - 1])) ||
                        !std::isdigit(static_cast<unsigned char>(value[j + 1])))
                        throw std::runtime_error("Invalid TOML number for '" + key + "'");
                }
                else number += c;
            }
            std::size_t consumed = 0;
            try
            {
                if (is_integer(key)) std::stoull(number, &consumed);
                else std::stod(number, &consumed);
            }
            catch (std::exception const &)
            {
                throw std::runtime_error("Invalid TOML number for '" + key + "'");
            }
            if (consumed != number.size())
                throw std::runtime_error("Invalid TOML number for '" + key + "'");
            value = number;
        }
        std::replace(key.begin(), key.end(), '_', '-');
        args.push_back("--" + key);
        args.push_back(value);
    }
    if (!input.eof()) throw std::runtime_error("Failed reading TOML config " + path.string());
    if (section != "v2") throw std::runtime_error("Missing [v2] TOML section in " + path.string());
    return args;
}

Options parse_options(int argc, char ** argv)
{
    Options o;
    std::vector<std::string> args{argv[0]};
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "--config" || (argc == 2 && arg.ends_with(".toml")))
        {
            if (!o.config_path.empty()) throw std::runtime_error("Only one --config is allowed");
            if (arg == "--config")
            {
                if (++i == argc) throw std::runtime_error("Missing value for --config");
                o.config_path = argv[i];
            }
            else o.config_path = arg;
        }
        else args.push_back(arg);
    }
    if (std::find(args.begin() + 1, args.end(), "--help") != args.end())
    {
        o.config_path.clear();
        args = {argv[0], "--help"};
    }
    else if (!o.config_path.empty())
    {
        auto config_args = load_v2_config_args(o.config_path);
        args.insert(args.begin() + 1, config_args.begin(), config_args.end());
    }
    for (std::size_t i = 1; i < args.size(); ++i)
    {
        std::string key = args[i];
        if (key == "--help")
        {
            std::cout << "Usage: lncrna_mers_v2 [--config FILE.toml] --reference REF.fa[.gz] "
                         "--queries Q.fa[.gz] --output OUT.tsv "
                         "[--k 15] [--mode ibf|exact|validate|calibrate] [--strand canonical|forward|both] "
                         "[--ambiguity skip|error] [--repeat include|exclude|annotate] "
                         "[--reference-repeat-bed REF.bed] [--query-repeat-bed QUERY.bed] "
                         "[--mask-target query|reference|both] [--repeat-classes OUT.tsv] "
                         "[--reference-hit-classes OUT.tsv] [--fragments OUT.tsv] "
                         "[--max-repeat-intervals 5000000] "
                         "[--exon-blocks BLOCKS.tsv] [--max-exon-blocks 1000000] "
                         "[--discover-loci OUT.tsv --locus-gap 8000 --max-locus-span 1000000] "
                         "[--max-locus-rows 100000] [--expected-loci-bed BED4/5] "
                         "[--expected-loci-output OUT.tsv] [--max-expected-loci 1000000] "
                         "[--max-query-bases 50000000] "
                         "[--fragment-size 8000] [--fpr 0.01] [--genome-fpr 0.01] [--hash-functions 3] "
                         "[--batch-windows 100000] [--basic-threshold 0.75] [--strict-threshold 0.85] "
                         "[--threshold-table TABLE.tsv --evaluation-control-set LABEL] "
                         "[--local-window BASES] [--positions OUT.tsv[.gz]] [--max-position-rows 1000000] "
                         "[--calibration-bins OUT.tsv] [--max-calibration-pairs 10000000] "
                         "[--max-ibf-bin-pairs 10000000] "
                         "[--save-index INDEX.ibf | --load-index INDEX.ibf] "
                         "[--metrics OUT.json]\n";
            std::exit(0);
        }
        if (i + 1 >= args.size()) throw std::runtime_error("Missing value for " + key);
        std::string value = args[++i];
        if (key == "--reference") o.reference = value;
        else if (key == "--queries") o.queries = value;
        else if (key == "--output") o.output = value;
        else if (key == "--metrics") o.metrics = value;
        else if (key == "--positions") o.positions = value;
        else if (key == "--calibration-bins") o.calibration_bins = value;
        else if (key == "--fragments") o.fragments = value;
        else if (key == "--reference-repeat-bed") o.reference_repeat_bed = value;
        else if (key == "--query-repeat-bed") o.query_repeat_bed = value;
        else if (key == "--repeat-classes") o.repeat_classes = value;
        else if (key == "--reference-hit-classes") o.reference_hit_classes = value;
        else if (key == "--exon-blocks") o.exon_blocks = value;
        else if (key == "--discover-loci") o.loci_output = value;
        else if (key == "--expected-loci-bed") o.expected_loci_bed = value;
        else if (key == "--expected-loci-output") o.expected_loci_output = value;
        else if (key == "--save-index") o.save_index = value;
        else if (key == "--load-index") o.load_index = value;
        else if (key == "--mask-target") o.mask_target = value;
        else if (key == "--threshold-table") o.threshold_table = value;
        else if (key == "--evaluation-control-set") o.evaluation_control_set = value;
        else if (key == "--k") o.k = std::stoull(value);
        else if (key == "--fragment-size") o.fragment_size = std::stoull(value);
        else if (key == "--hash-functions") o.hash_functions = std::stoull(value);
        else if (key == "--batch-windows") o.batch_windows = std::stoull(value);
        else if (key == "--local-window") { o.local_window = std::stoull(value); o.local_window_given = true; }
        else if (key == "--max-position-rows") o.max_position_rows = std::stoull(value);
        else if (key == "--max-calibration-pairs") o.max_calibration_pairs = std::stoull(value);
        else if (key == "--max-ibf-bin-pairs") o.max_ibf_bin_pairs = std::stoull(value);
        else if (key == "--max-repeat-intervals") o.max_repeat_intervals = std::stoull(value);
        else if (key == "--max-exon-blocks") o.max_exon_blocks = std::stoull(value);
        else if (key == "--max-query-bases") o.max_query_bases = std::stoull(value);
        else if (key == "--locus-gap") o.locus_gap = std::stoull(value);
        else if (key == "--max-locus-span") o.max_locus_span = std::stoull(value);
        else if (key == "--max-locus-rows") o.max_locus_rows = std::stoull(value);
        else if (key == "--max-expected-loci") o.max_expected_loci = std::stoull(value);
        else if (key == "--fpr") o.fpr = std::stod(value);
        else if (key == "--genome-fpr") { o.genome_fpr = std::stod(value); o.genome_fpr_given = true; }
        else if (key == "--basic-threshold") o.basic = std::stod(value);
        else if (key == "--strict-threshold") o.strict = std::stod(value);
        else if (key == "--strand") o.strand = value;
        else if (key == "--ambiguity") o.ambiguity = value;
        else if (key == "--repeat") o.repeat = value;
        else if (key == "--mode") o.mode = value;
        else throw std::runtime_error("Unknown option: " + key);
    }
    if (o.reference.empty() || o.queries.empty() || o.output.empty())
        throw std::runtime_error("--reference, --queries, and --output are required");
    if (o.k < 1 || o.k > 32 || o.fragment_size < o.k || o.hash_functions < 1 ||
        o.hash_functions > 32 || o.batch_windows < 1 || o.max_position_rows < 1 ||
        o.max_calibration_pairs < 1 || o.max_ibf_bin_pairs < 1 ||
        o.max_repeat_intervals < 1 || o.max_exon_blocks < 1 || o.max_query_bases < o.k ||
        o.max_locus_span < o.k || o.max_locus_rows < 1 || o.max_expected_loci < 1 ||
        !(o.fpr > 0 && o.fpr <= 0.5) || !(o.genome_fpr > 0 && o.genome_fpr <= 0.5) ||
        !(0 <= o.basic && o.basic <= o.strict && o.strict <= 1))
        throw std::runtime_error("Invalid numeric option (k 1..32, fragment >= k, FPR (0,0.5], thresholds 0..1)");
    if (o.strand != "canonical" && o.strand != "forward" && o.strand != "both")
        throw std::runtime_error("Invalid --strand");
    if (o.ambiguity != "skip" && o.ambiguity != "error") throw std::runtime_error("Invalid --ambiguity");
    if (o.repeat != "include" && o.repeat != "exclude" && o.repeat != "annotate")
        throw std::runtime_error("Invalid --repeat");
    if (o.mask_target != "query" && o.mask_target != "reference" && o.mask_target != "both")
        throw std::runtime_error("Invalid --mask-target");
    if (o.repeat != "annotate" && !o.repeat_classes.empty())
        throw std::runtime_error("--repeat-classes requires --repeat annotate");
    if (!o.reference_hit_classes.empty() && o.repeat != "annotate")
        throw std::runtime_error("--reference-hit-classes requires --repeat annotate");
    if (o.mode != "ibf" && o.mode != "exact" && o.mode != "validate" && o.mode != "calibrate")
        throw std::runtime_error("Invalid --mode");
    if (o.mode == "ibf" && (o.local_window || !o.positions.empty() || !o.loci_output.empty() ||
                           !o.expected_loci_bed.empty()))
        throw std::runtime_error("Localization and exact positions require --mode validate, calibrate, or exact");
    if (o.expected_loci_bed.empty() && !o.expected_loci_output.empty())
        throw std::runtime_error("--expected-loci-output requires --expected-loci-bed");
    if (o.mode == "ibf" && !o.exon_blocks.empty())
        throw std::runtime_error("--exon-blocks requires --mode validate, calibrate, or exact");
    if (o.mode == "ibf" && !o.reference_hit_classes.empty())
        throw std::runtime_error("Reference hit classes require exact positions");
    if (o.mode == "exact" && o.genome_fpr_given)
        throw std::runtime_error("--genome-fpr requires an IBF mode");
    if (o.mode != "calibrate" && !o.calibration_bins.empty())
        throw std::runtime_error("--calibration-bins requires --mode calibrate");
    if (!o.save_index.empty() && !o.load_index.empty())
        throw std::runtime_error("Choose --save-index or --load-index, not both");
    if (o.mode == "exact" && (!o.save_index.empty() || !o.load_index.empty()))
        throw std::runtime_error("Index serialization requires an IBF mode");
    if (o.local_window_given && o.local_window < o.k)
        throw std::runtime_error("--local-window must be at least k");
    if (o.evaluation_control_set.find_first_of("\t\r\n") != std::string::npos)
        throw std::runtime_error("Invalid evaluation control-set label");
    load_thresholds(o);
    o.requested_per_bin_fpr = o.fpr;
    if (o.metrics.empty()) o.metrics = o.output.string() + ".metrics.json";
    if (o.mode == "calibrate" && o.calibration_bins.empty()) o.calibration_bins = o.output.string() + ".bins.tsv";
    if (o.repeat == "annotate" && o.repeat_classes.empty()) o.repeat_classes = o.output.string() + ".repeat_classes.tsv";
    if (!o.expected_loci_bed.empty() && o.expected_loci_output.empty())
        o.expected_loci_output = o.output.string() + ".expected_loci.tsv";
    std::vector<fs::path> inputs{o.reference, o.queries}, outputs{o.output, o.metrics};
    if (!o.config_path.empty()) inputs.push_back(o.config_path);
    for (auto const & input : {o.reference_repeat_bed, o.query_repeat_bed, o.exon_blocks, o.expected_loci_bed,
                              o.threshold_table,
                              o.load_index})
        if (!input.empty()) inputs.push_back(input);
    for (auto const & output : {o.positions, o.calibration_bins, o.fragments, o.repeat_classes,
                               o.reference_hit_classes, o.loci_output,
                               o.expected_loci_output, o.save_index})
        if (!output.empty()) outputs.push_back(output);
    if (!o.save_index.empty()) outputs.push_back(o.save_index.string() + ".metadata.json");
    auto same_path = [](fs::path const & a, fs::path const & b)
    {
        return fs::absolute(a).lexically_normal() == fs::absolute(b).lexically_normal() ||
               (fs::exists(a) && fs::exists(b) && fs::equivalent(a, b));
    };
    for (std::size_t i = 0; i < outputs.size(); ++i)
    {
        for (auto const & input : inputs)
            if (same_path(outputs[i], input)) throw std::runtime_error("Input and output paths must be distinct");
        for (std::size_t j = i + 1; j < outputs.size(); ++j)
            if (same_path(outputs[i], outputs[j])) throw std::runtime_error("Input and output paths must be distinct");
    }
    return o;
}

// Buffered byte reader keeps memory independent of FASTA line or contig length.
class ByteReader
{
public:
    explicit ByteReader(fs::path const & path) : gzip_{path.extension() == ".gz"}
    {
        if (gzip_)
        {
            gz_ = gzopen(path.string().c_str(), "rb");
            if (!gz_) throw std::runtime_error("Cannot open " + path.string());
        }
        else
        {
            file_.open(path, std::ios::binary);
            if (!file_) throw std::runtime_error("Cannot open " + path.string());
        }
    }
    ~ByteReader() { if (gz_) gzclose(gz_); }
    int get()
    {
        if (at_ == length_)
        {
            at_ = 0;
            if (gzip_)
            {
                int n = gzread(gz_, buffer_, sizeof(buffer_));
                if (n < 0) throw std::runtime_error("gzip read failed");
                length_ = static_cast<std::size_t>(n);
            }
            else
            {
                file_.read(buffer_, sizeof(buffer_));
                length_ = static_cast<std::size_t>(file_.gcount());
            }
            if (!length_) return -1;
        }
        return static_cast<unsigned char>(buffer_[at_++]);
    }
private:
    bool gzip_;
    gzFile gz_ = nullptr;
    std::ifstream file_;
    char buffer_[65536]{};
    std::size_t at_ = 0, length_ = 0;
};

template <typename OnStart, typename OnBase>
void scan_fasta(fs::path const & path, OnStart on_start, OnBase on_base)
{
    ByteReader in(path);
    bool line_start = true, in_header = false, seen_header = false;
    std::string header;
    std::size_t pos = 0;
    int raw;
    while ((raw = in.get()) != -1)
    {
        char c = static_cast<char>(raw);
        if (c == '\n' || c == '\r')
        {
            if (in_header)
            {
                std::string id = header.substr(0, header.find_first_of(" \t"));
                if (id.empty()) throw std::runtime_error("Empty FASTA identifier in " + path.string());
                on_start(id);
                seen_header = true;
                pos = 0;
                header.clear();
                in_header = false;
            }
            line_start = true;
            continue;
        }
        if (line_start && c == '>')
        {
            in_header = true;
            header.clear();
            line_start = false;
            continue;
        }
        line_start = false;
        if (in_header) header += c;
        else if (c != ' ' && c != '\t')
        {
            if (!seen_header) throw std::runtime_error("Sequence before FASTA header in " + path.string());
            on_base(pos++, c);
        }
    }
    if (in_header)
    {
        std::string id = header.substr(0, header.find_first_of(" \t"));
        if (id.empty()) throw std::runtime_error("Empty FASTA identifier in " + path.string());
        on_start(id);
        seen_header = true;
    }
    if (!seen_header) throw std::runtime_error("No FASTA records in " + path.string());
}

constexpr uint32_t mixed_repeat_class = std::numeric_limits<uint32_t>::max();

struct RepeatInterval { std::size_t start, end; uint32_t class_id; };
struct RepeatBed
{
    std::unordered_map<std::string, std::vector<RepeatInterval>> intervals;
    std::vector<std::string> classes{"NONE", "SOFTMASK"};
};

RepeatBed load_repeat_bed(fs::path const & path, std::size_t max_intervals)
{
    RepeatBed bed;
    if (path.empty()) return bed;
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open repeat BED " + path.string());
    std::unordered_map<std::string, uint32_t> class_ids;
    std::string line;
    std::size_t count = 0, line_no = 0;
    while (std::getline(input, line))
    {
        ++line_no;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        auto fields = split_tabs(line);
        if (fields.size() < 3 || fields[0].empty())
            throw std::runtime_error("Invalid repeat BED row " + std::to_string(line_no));
        std::size_t start, end;
        try
        {
            std::size_t n1 = 0, n2 = 0;
            if (fields[1].empty() || fields[2].empty() || fields[1][0] == '-' || fields[2][0] == '-')
                throw std::invalid_argument{"invalid interval"};
            start = std::stoull(fields[1], &n1);
            end = std::stoull(fields[2], &n2);
            if (n1 != fields[1].size() || n2 != fields[2].size() || start >= end)
                throw std::invalid_argument{"invalid interval"};
        }
        catch (std::exception const &)
        {
            throw std::runtime_error("Invalid repeat BED coordinates on row " + std::to_string(line_no));
        }
        std::string name = fields.size() > 3 && !fields[3].empty() ? fields[3] : "UNCLASSIFIED";
        if (name == "NONE" || name == "SOFTMASK" || name == "MULTIPLE_CLASSES")
            throw std::runtime_error("Reserved repeat class on row " + std::to_string(line_no));
        auto [it, inserted] = class_ids.try_emplace(name, static_cast<uint32_t>(bed.classes.size()));
        if (inserted) bed.classes.push_back(name);
        bed.intervals[fields[0]].push_back({start, end, it->second});
        if (++count > max_intervals) throw std::runtime_error("Repeat BED exceeds --max-repeat-intervals");
    }
    for (auto & [contig, intervals] : bed.intervals)
        std::sort(intervals.begin(), intervals.end(), [](auto const & a, auto const & b)
        { return std::tie(a.start, a.end, a.class_id) < std::tie(b.start, b.end, b.class_id); });
    return bed;
}

class BedCursor
{
public:
    explicit BedCursor(RepeatBed const & bed) : bed_{bed} {}
    void reset(std::string const & id)
    {
        auto it = bed_.intervals.find(id);
        intervals_ = it == bed_.intervals.end() ? nullptr : &it->second;
        next_ = 0;
        active_.clear();
    }
    uint32_t at(std::size_t pos)
    {
        if (!intervals_) return 0;
        while (next_ < intervals_->size() && (*intervals_)[next_].start <= pos)
        {
            if ((*intervals_)[next_].end > pos) active_.push_back(next_);
            ++next_;
        }
        active_.erase(std::remove_if(active_.begin(), active_.end(), [&](std::size_t i)
        { return (*intervals_)[i].end <= pos; }), active_.end());
        uint32_t id = 0;
        for (std::size_t i : active_)
        {
            uint32_t current = (*intervals_)[i].class_id;
            if (id && id != current) return mixed_repeat_class;
            id = current;
        }
        return id;
    }
private:
    RepeatBed const & bed_;
    std::vector<RepeatInterval> const * intervals_ = nullptr;
    std::vector<std::size_t> active_;
    std::size_t next_ = 0;
};

struct Window
{
    uint64_t forward = 0, reverse = 0;
    bool valid = false, repeat = false, base_valid = false;
    uint32_t repeat_class = 0;
};

class WindowMaker
{
public:
    explicit WindowMaker(Options const & o, bool exclude_here) : o_{o}, exclude_here_{exclude_here}, core_{o.k} {}
    void reset() { codes_.clear(); masked_ = 0; core_.reset(); }
    std::optional<Window> push(char c, uint32_t bed_class = 0)
    {
        if (fraglncrna::dna_code(c) < 0 && o_.ambiguity == "error")
            throw std::runtime_error(std::string("Ambiguous base in FASTA: ") + c);
        auto core = core_.push(c);
        bool softmask = c >= 'a' && c <= 'z';
        bool repeat = softmask || bed_class != 0;
        if (codes_.size() == o_.k)
        {
            masked_ -= codes_.front().repeat;
            codes_.pop_front();
        }
        codes_.push_back({repeat, softmask ? uint32_t{1} : uint32_t{0}, bed_class});
        masked_ += repeat;
        if (!core) return std::nullopt;
        bool base_valid = core->valid;
        uint32_t class_id = 0;
        if (masked_)
            for (auto const & tag : codes_)
                for (uint32_t current : {tag.soft_class, tag.bed_class})
                    if (current)
                        class_id = class_id && class_id != current ? mixed_repeat_class : current;
        return Window{core->forward, core->reverse, base_valid && (!exclude_here_ || masked_ == 0),
                      masked_ != 0, base_valid, class_id};
    }
private:
    struct BaseTag { bool repeat; uint32_t soft_class, bed_class; };
    Options const & o_;
    bool exclude_here_;
    fraglncrna::RollingDnaWindow core_;
    std::deque<BaseTag> codes_;
    std::size_t masked_ = 0;
};

uint64_t index_key(Window const & w, Options const & o)
{
    return o.strand == "canonical" ? std::min(w.forward, w.reverse) : w.forward;
}

struct ReferenceInfo
{
    std::vector<std::string> ids;
    std::vector<std::size_t> lengths, bin_offsets;
    std::vector<std::size_t> bin_insertions, bin_distinct;
    std::size_t total_bases = 0, valid_windows = 0, bin_count = 0;
    std::size_t bin_bits = 0;
    uint64_t sequence_hash = 14695981039346656037ULL;
    std::string sequence_sha256;
};

struct ExonBlock
{
    std::size_t query_start, query_end, contig, ref_start, ref_end;
    char strand;
};

using ExonBlocks = std::unordered_map<std::string, std::vector<ExonBlock>>;

struct ExpectedLocus
{
    std::size_t contig, start, end;
    std::string label;
};
using ExpectedLoci = std::unordered_map<std::string, std::vector<ExpectedLocus>>;

std::size_t parse_coordinate(std::string const & value, std::string const & context)
{
    if (value.empty() || value[0] == '-') throw std::runtime_error("Invalid coordinate in " + context);
    try
    {
        std::size_t parsed = 0;
        auto result = std::stoull(value, &parsed);
        if (parsed == value.size()) return result;
    }
    catch (std::exception const &) {}
    throw std::runtime_error("Invalid coordinate in " + context);
}

ExonBlocks load_exon_blocks(Options const & o, ReferenceInfo const & info)
{
    ExonBlocks blocks;
    if (o.exon_blocks.empty()) return blocks;
    std::ifstream input(o.exon_blocks);
    if (!input) throw std::runtime_error("Cannot open exon blocks " + o.exon_blocks.string());
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("Empty exon-block TSV");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line != "query_id\tquery_start\tquery_end\tcontig\tref_start\tref_end\tstrand")
        throw std::runtime_error("Invalid exon-block TSV header");
    std::unordered_map<std::string, std::size_t> contigs;
    for (std::size_t i = 0; i < info.ids.size(); ++i) contigs[info.ids[i]] = i;
    std::size_t count = 0, line_no = 1;
    while (std::getline(input, line))
    {
        ++line_no;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        auto fields = split_tabs(line);
        if (fields.size() != 7 || fields[0].empty() || !contigs.contains(fields[3]) ||
            (fields[6] != "+" && fields[6] != "-"))
            throw std::runtime_error("Invalid exon-block row " + std::to_string(line_no));
        std::string context = "exon-block row " + std::to_string(line_no);
        std::size_t qs = parse_coordinate(fields[1], context), qe = parse_coordinate(fields[2], context);
        std::size_t rs = parse_coordinate(fields[4], context), re = parse_coordinate(fields[5], context);
        if (qs >= qe || rs >= re || qe - qs != re - rs || re > info.lengths[contigs.at(fields[3])])
            throw std::runtime_error("Invalid exon-block span on row " + std::to_string(line_no));
        blocks[fields[0]].push_back({qs, qe, contigs.at(fields[3]), rs, re, fields[6][0]});
        if (++count > o.max_exon_blocks) throw std::runtime_error("Exon blocks exceed --max-exon-blocks");
    }
    for (auto & [id, query_blocks] : blocks)
    {
        std::sort(query_blocks.begin(), query_blocks.end(), [](auto const & a, auto const & b)
        { return a.query_start < b.query_start; });
        for (std::size_t i = 1; i < query_blocks.size(); ++i)
            if (query_blocks[i - 1].query_end > query_blocks[i].query_start)
                throw std::runtime_error("Overlapping query exon blocks for " + id);
    }
    return blocks;
}

ExpectedLoci load_expected_loci(Options const & o, ReferenceInfo const & info)
{
    ExpectedLoci loci;
    if (o.expected_loci_bed.empty()) return loci;
    std::ifstream input(o.expected_loci_bed);
    if (!input) throw std::runtime_error("Cannot open expected-locus BED " + o.expected_loci_bed.string());
    std::unordered_map<std::string, std::size_t> contigs;
    for (std::size_t i = 0; i < info.ids.size(); ++i) contigs[info.ids[i]] = i;
    std::string line;
    std::size_t count = 0, line_no = 0;
    while (std::getline(input, line))
    {
        ++line_no;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        auto fields = split_tabs(line);
        if ((fields.size() != 4 && fields.size() != 5) || !contigs.contains(fields[0]) || fields[3].empty())
            throw std::runtime_error("Invalid expected-locus BED row " + std::to_string(line_no));
        std::string context = "expected-locus BED row " + std::to_string(line_no);
        std::size_t start = parse_coordinate(fields[1], context), end = parse_coordinate(fields[2], context);
        std::size_t contig = contigs.at(fields[0]);
        if (start >= end || end > info.lengths[contig])
            throw std::runtime_error("Invalid expected-locus interval on row " + std::to_string(line_no));
        std::string label = fields.size() == 5 && !fields[4].empty() ? fields[4] : ".";
        loci[fields[3]].push_back({contig, start, end, label});
        if (++count > o.max_expected_loci) throw std::runtime_error("Expected loci exceed --max-expected-loci");
    }
    for (auto & [id, regions] : loci)
    {
        std::sort(regions.begin(), regions.end(), [](auto const & a, auto const & b)
        { return std::tie(a.contig, a.start, a.end) < std::tie(b.contig, b.start, b.end); });
        for (std::size_t i = 1; i < regions.size(); ++i)
            if (regions[i - 1].contig == regions[i].contig && regions[i - 1].end > regions[i].start)
                throw std::runtime_error("Overlapping expected loci for query " + id);
    }
    return loci;
}

ReferenceInfo inspect_reference(Options const & o, RepeatBed const & reference_bed)
{
    ReferenceInfo info;
    fraglncrna::Sha256 sha;
    sha.update("FragIncRNA-FASTA-semantic-v1");
    WindowMaker maker{o, o.repeat == "exclude" && o.mask_target != "query"};
    BedCursor bed{reference_bed};
    std::unordered_set<std::string> seen_ids;
    scan_fasta(o.reference,
        [&](std::string const & id)
        {
            if (!seen_ids.insert(id).second)
                throw std::runtime_error("Duplicate reference FASTA identifier: " + id);
            if (!info.ids.empty()) sha256_record_end(sha, info.lengths.back());
            info.sequence_hash = fnv_update(info.sequence_hash, 255);
            sha256_record_start(sha, id);
            for (unsigned char c : id) info.sequence_hash = fnv_update(info.sequence_hash, c);
            info.ids.push_back(id);
            info.bin_offsets.push_back(info.bin_count);
            info.lengths.push_back(0);
            maker.reset();
            bed.reset(id);
        },
        [&](std::size_t pos, char c)
        {
            info.sequence_hash = fnv_update(info.sequence_hash, static_cast<unsigned char>(c));
            sha.update(static_cast<unsigned char>(c));
            ++info.lengths.back();
            ++info.total_bases;
            if (auto w = maker.push(c, bed.at(pos)); w && w->valid) ++info.valid_windows;
        });
    for (std::size_t i = 0; i < info.lengths.size(); ++i)
    {
        info.bin_offsets[i] = info.bin_count;
        std::size_t bins = (info.lengths[i] + o.fragment_size - 1) / o.fragment_size;
        info.bin_count += bins;
    }
    info.bin_insertions.resize(info.bin_count);
    info.bin_distinct.resize(info.bin_count);
    sha256_record_end(sha, info.lengths.back());
    info.sequence_sha256 = sha.hex_digest();
    return info;
}

using Ibf = seqan3::interleaved_bloom_filter<>;

std::unique_ptr<Ibf> build_index(Options const & o, ReferenceInfo & info, RepeatBed const & reference_bed)
{
    double n = static_cast<double>(o.fragment_size);
    double h = static_cast<double>(o.hash_functions);
    std::size_t bits = std::max<std::size_t>(1024, static_cast<std::size_t>(
        std::ceil(-h * n / std::log(1.0 - std::pow(o.fpr, 1.0 / h)))));
    info.bin_bits = bits;
    auto ibf = std::make_unique<Ibf>(seqan3::bin_count{std::max<std::size_t>(1, info.bin_count)},
                                     seqan3::bin_size{bits}, seqan3::hash_function_count{o.hash_functions});
    WindowMaker maker{o, o.repeat == "exclude" && o.mask_target != "query"};
    BedCursor bed{reference_bed};
    std::size_t record = 0;
    std::size_t active_bin = std::numeric_limits<std::size_t>::max();
    std::unordered_set<uint64_t> distinct;
    scan_fasta(o.reference,
        [&](std::string const & id) { maker.reset(); bed.reset(id); ++record; },
        [&](std::size_t pos, char c)
        {
            auto w = maker.push(c, bed.at(pos));
            if (!w || !w->valid) return;
            std::size_t start = pos + 1 - o.k;
            std::size_t bin = info.bin_offsets[record - 1] + start / o.fragment_size;
            if (bin != active_bin)
            {
                if (active_bin != std::numeric_limits<std::size_t>::max())
                    info.bin_distinct[active_bin] = distinct.size();
                distinct.clear();
                active_bin = bin;
            }
            uint64_t key = index_key(*w, o);
            distinct.insert(key);
            ibf->emplace(key, seqan3::bin_index{bin});
            ++info.bin_insertions[bin];
        });
    if (active_bin != std::numeric_limits<std::size_t>::max())
        info.bin_distinct[active_bin] = distinct.size();
    return ibf;
}

struct IndexMetadata
{
    uint32_t format_version = 3;
    std::string software_version = "2-preview.3";
    std::string compiler_version = __VERSION__;
    std::string seqan_version = "3.4.1";
    uint32_t size_t_bytes = sizeof(std::size_t);
    std::string byte_order = std::endian::native == std::endian::little ? "little" : "big";
    std::string reference_sha256, reference_repeat_bed_sha256;
    std::string strand, ambiguity, repeat, mask_target;
    uint64_t k = 0, fragment_size = 0, hash_functions = 0, bin_count = 0, bin_bits = 0;
    uint64_t valid_reference_windows = 0;
    double configured_fpr = 0;
    std::vector<std::string> contig_ids;
    std::vector<std::size_t> contig_lengths, bin_insertions, bin_distinct;

    template <class Archive>
    void serialize(Archive & archive)
    {
        archive(format_version, software_version, compiler_version, seqan_version, size_t_bytes, byte_order,
                reference_sha256, reference_repeat_bed_sha256, strand, ambiguity, repeat, mask_target,
                k, fragment_size, hash_functions, bin_count, bin_bits, valid_reference_windows, configured_fpr,
                contig_ids, contig_lengths, bin_insertions, bin_distinct);
    }
};

IndexMetadata make_index_metadata(Options const & o, ReferenceInfo const & info,
                                  std::string const & reference_bed_sha)
{
    IndexMetadata meta;
    meta.reference_sha256 = info.sequence_sha256;
    meta.reference_repeat_bed_sha256 = reference_bed_sha;
    meta.strand = o.strand;
    meta.ambiguity = o.ambiguity;
    meta.repeat = o.repeat;
    meta.mask_target = o.mask_target;
    meta.k = o.k;
    meta.fragment_size = o.fragment_size;
    meta.hash_functions = o.hash_functions;
    meta.bin_count = info.bin_count;
    meta.bin_bits = info.bin_bits;
    meta.valid_reference_windows = info.valid_windows;
    meta.configured_fpr = o.fpr;
    meta.contig_ids = info.ids;
    meta.contig_lengths = info.lengths;
    meta.bin_insertions = info.bin_insertions;
    meta.bin_distinct = info.bin_distinct;
    return meta;
}

void check_index_metadata(IndexMetadata const & stored, IndexMetadata const & expected)
{
    auto require_same = [](auto const & a, auto const & b, char const * field)
    {
        if (a != b) throw std::runtime_error(std::string("Incompatible index metadata: ") + field);
    };
    require_same(stored.format_version, expected.format_version, "format_version");
    require_same(stored.software_version, expected.software_version, "software_version");
    require_same(stored.compiler_version, expected.compiler_version, "compiler_version");
    require_same(stored.seqan_version, expected.seqan_version, "seqan_version");
    require_same(stored.size_t_bytes, expected.size_t_bytes, "size_t_bytes");
    require_same(stored.byte_order, expected.byte_order, "byte_order");
    require_same(stored.reference_sha256, expected.reference_sha256, "reference_sha256");
    require_same(stored.reference_repeat_bed_sha256, expected.reference_repeat_bed_sha256,
                 "reference_repeat_bed_sha256");
    require_same(stored.strand, expected.strand, "strand");
    require_same(stored.ambiguity, expected.ambiguity, "ambiguity");
    require_same(stored.repeat, expected.repeat, "repeat");
    require_same(stored.mask_target, expected.mask_target, "mask_target");
    require_same(stored.k, expected.k, "k");
    require_same(stored.fragment_size, expected.fragment_size, "fragment_size");
    require_same(stored.hash_functions, expected.hash_functions, "hash_functions");
    require_same(stored.bin_count, expected.bin_count, "bin_count");
    require_same(stored.valid_reference_windows, expected.valid_reference_windows, "valid_reference_windows");
    require_same(stored.configured_fpr, expected.configured_fpr, "configured_fpr");
    require_same(stored.contig_ids, expected.contig_ids, "contig_ids");
    require_same(stored.contig_lengths, expected.contig_lengths, "contig_lengths");
    if (stored.bin_insertions.size() != stored.bin_count || stored.bin_distinct.size() != stored.bin_count ||
        stored.bin_bits == 0)
        throw std::runtime_error("Incompatible index metadata: bin layout");
    for (std::size_t bin = 0; bin < stored.bin_count; ++bin)
        if (stored.bin_distinct[bin] > stored.bin_insertions[bin])
            throw std::runtime_error("Incompatible index metadata: distinct-key count");
}

constexpr std::string_view index_magic = "FRAGLNCRNA-IBF3";

void save_index_file(fs::path const & path, IndexMetadata const & meta, Ibf const & ibf)
{
    fs::create_directories(path.parent_path().empty() ? fs::path{"."} : path.parent_path());
    std::ostringstream header_stream(std::ios::out | std::ios::binary);
    { cereal::BinaryOutputArchive archive(header_stream); archive(meta); }
    std::string header = header_stream.str();
    if (header.size() > 64 * 1024 * 1024) throw std::runtime_error("Index metadata exceeds 64 MiB");
    fs::path temporary = path.string() + ".partial." + std::to_string(getpid());
    try
    {
        std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("Cannot write index " + temporary.string());
        output.write(index_magic.data(), static_cast<std::streamsize>(index_magic.size()));
        uint64_t header_size = header.size();
        output.write(reinterpret_cast<char const *>(&header_size), sizeof(header_size));
        std::string header_sha = sha256_text(header);
        output.write(header_sha.data(), static_cast<std::streamsize>(header_sha.size()));
        output.write(header.data(), static_cast<std::streamsize>(header.size()));
        std::uint64_t payload_digest_offset = index_magic.size() + sizeof(header_size) +
                                              header_sha.size() + header.size();
        std::string placeholder(64, '0');
        output.write(placeholder.data(), static_cast<std::streamsize>(placeholder.size()));
        { cereal::BinaryOutputArchive archive(output); archive(ibf); }
        output.close();
        if (!output) throw std::runtime_error("Failed writing index " + temporary.string());
        std::string payload_sha = sha256_file(temporary, payload_digest_offset + placeholder.size());
        std::fstream patch(temporary, std::ios::binary | std::ios::in | std::ios::out);
        if (!patch) throw std::runtime_error("Cannot finalize index " + temporary.string());
        patch.seekp(static_cast<std::streamoff>(payload_digest_offset));
        patch.write(payload_sha.data(), static_cast<std::streamsize>(payload_sha.size()));
        patch.close();
        if (!patch) throw std::runtime_error("Failed finalizing index " + temporary.string());
        fs::rename(temporary, path);
    }
    catch (...)
    {
        std::error_code ignored;
        fs::remove(temporary, ignored);
        throw;
    }
}

std::unique_ptr<Ibf> load_index_file(fs::path const & path, IndexMetadata const & expected,
                                     ReferenceInfo & info)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open index " + path.string());
    std::string magic(index_magic.size(), '\0');
    input.read(magic.data(), static_cast<std::streamsize>(magic.size()));
    if (magic != index_magic) throw std::runtime_error("Incompatible index format: bad magic");
    uint64_t header_size = 0;
    input.read(reinterpret_cast<char *>(&header_size), sizeof(header_size));
    if (!input || header_size == 0 || header_size > 64 * 1024 * 1024 ||
        header_size > fs::file_size(path))
        throw std::runtime_error("Incompatible index format: invalid metadata size");
    std::string expected_header_sha(64, '\0');
    input.read(expected_header_sha.data(), 64);
    if (!input) throw std::runtime_error("Truncated index metadata checksum");
    std::string header(static_cast<std::size_t>(header_size), '\0');
    input.read(header.data(), static_cast<std::streamsize>(header.size()));
    if (!input) throw std::runtime_error("Truncated index metadata");
    if (sha256_text(header) != expected_header_sha)
        throw std::runtime_error("Index metadata checksum mismatch");
    std::string expected_payload_sha(64, '\0');
    input.read(expected_payload_sha.data(), 64);
    if (!input) throw std::runtime_error("Truncated index payload checksum");
    std::uint64_t payload_offset = index_magic.size() + sizeof(header_size) + 64 + header_size + 64;
    if (sha256_file(path, payload_offset) != expected_payload_sha)
        throw std::runtime_error("Index payload checksum mismatch");
    std::istringstream header_input(header, std::ios::in | std::ios::binary);
    IndexMetadata stored;
    { cereal::BinaryInputArchive archive(header_input); archive(stored); }
    check_index_metadata(stored, expected);
    auto ibf = std::make_unique<Ibf>();
    { cereal::BinaryInputArchive archive(input); archive(*ibf); }
    if (!input || ibf->bin_count() != std::max<std::size_t>(1, info.bin_count) ||
        ibf->hash_function_count() != expected.hash_functions || ibf->bin_size() != stored.bin_bits)
        throw std::runtime_error("Incompatible or truncated index payload");
    info.bin_bits = stored.bin_bits;
    info.bin_insertions = std::move(stored.bin_insertions);
    info.bin_distinct = std::move(stored.bin_distinct);
    return ibf;
}

void write_index_metadata_json(fs::path const & index_path, IndexMetadata const & meta,
                               std::string const & index_sha, std::size_t disk_bytes)
{
    fs::path path = index_path.string() + ".metadata.json";
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write index metadata " + path.string());
    out << "{\n  \"index_format_version\": " << meta.format_version
        << ",\n  \"semantic_fasta_hash_encoding\": \"FragIncRNA-FASTA-semantic-v1\""
        << ",\n  \"software_version\": " << json_string(meta.software_version)
        << ",\n  \"compiler_version\": " << json_string(meta.compiler_version)
        << ",\n  \"seqan_version\": " << json_string(meta.seqan_version)
        << ",\n  \"size_t_bytes\": " << meta.size_t_bytes
        << ",\n  \"byte_order\": " << json_string(meta.byte_order)
        << ",\n  \"reference_sequence_sha256\": " << json_string(meta.reference_sha256)
        << ",\n  \"reference_repeat_bed_sha256\": " << json_string(meta.reference_repeat_bed_sha256)
        << ",\n  \"strand_mode\": " << json_string(meta.strand)
        << ",\n  \"ambiguity_policy\": " << json_string(meta.ambiguity)
        << ",\n  \"repeat_policy\": " << json_string(meta.repeat)
        << ",\n  \"mask_target\": " << json_string(meta.mask_target)
        << ",\n  \"k\": " << meta.k << ",\n  \"fragment_size\": " << meta.fragment_size
        << ",\n  \"fragment_overlap\": " << (meta.k - 1)
        << ",\n  \"bin_count\": " << meta.bin_count
        << ",\n  \"bin_bits\": " << meta.bin_bits
        << ",\n  \"hash_functions\": " << meta.hash_functions
        << ",\n  \"configured_per_bin_fpr\": " << std::setprecision(17) << meta.configured_fpr
        << ",\n  \"valid_reference_windows\": " << meta.valid_reference_windows
        << ",\n  \"distinct_bin_kmers_sum\": "
        << std::accumulate(meta.bin_distinct.begin(), meta.bin_distinct.end(), std::size_t{0})
        << ",\n  \"disk_index_bytes\": " << disk_bytes
        << ",\n  \"index_file_sha256\": " << json_string(index_sha)
        << ",\n  \"contigs\": [";
    for (std::size_t i = 0; i < meta.contig_ids.size(); ++i)
    {
        if (i) out << ',';
        out << "\n    {\"id\": " << json_string(meta.contig_ids[i])
            << ", \"length\": " << meta.contig_lengths[i] << '}';
    }
    out << "\n  ]\n}\n";
    if (!out) throw std::runtime_error("Failed writing index metadata");
}

struct Query
{
    std::string id;
    std::size_t index = 0, total_windows = 0, valid_windows = 0, base_valid_windows = 0, repeat_windows = 0;
    std::size_t positive_bin_pairs = 0;
    std::vector<Window> windows;
    std::vector<std::size_t> bin_hits;
    std::vector<std::vector<std::size_t>> positive_bins;
};

struct Hit
{
    std::size_t window = 0, contig = 0, start = 0;
    char strand = '+';
    uint32_t reference_repeat_class = 0;
};

struct ExactResult
{
    std::unordered_set<uint64_t> found;
    std::vector<std::vector<Hit>> hits;
    std::unordered_map<uint64_t, std::vector<std::size_t>> exact_bins;
};

class PositionWriter
{
public:
    explicit PositionWriter(fs::path const & path) : compressed_{path.extension() == ".gz"}
    {
        if (path.empty()) return;
        fs::create_directories(path.parent_path().empty() ? fs::path{"."} : path.parent_path());
        if (compressed_)
        {
            gz_ = gzopen(path.string().c_str(), "wb");
            if (!gz_) throw std::runtime_error("Cannot write " + path.string());
        }
        else
        {
            file_.open(path);
            if (!file_) throw std::runtime_error("Cannot write " + path.string());
        }
        write("query_index\tquery_id\treference_id\tcontig\tstart\tend\tstrand\tquery_window_index\n");
    }
    ~PositionWriter() { if (gz_) gzclose(gz_); }
    void close()
    {
        if (gz_)
        {
            int result = gzclose(gz_);
            gz_ = nullptr;
            if (result != Z_OK) throw std::runtime_error("Failed closing compressed positions");
        }
        if (file_.is_open())
        {
            file_.close();
            if (!file_) throw std::runtime_error("Failed closing positions");
        }
    }
    void write(std::string const & line)
    {
        if (gz_)
        {
            if (gzwrite(gz_, line.data(), static_cast<unsigned>(line.size())) != static_cast<int>(line.size()))
                throw std::runtime_error("Failed writing compressed positions");
        }
        else if (file_)
        {
            file_ << line;
            if (!file_) throw std::runtime_error("Failed writing positions");
        }
    }
private:
    bool compressed_;
    gzFile gz_ = nullptr;
    std::ofstream file_;
};

std::size_t count_bins(Ibf::membership_agent_type & agent, Window const & w, Options const & o)
{
    auto const & first = agent.bulk_contains(index_key(w, o));
    if (o.strand != "both" || w.forward == w.reverse)
        return std::ranges::count(first, true);
    std::vector<bool> union_bins(first.begin(), first.end());
    auto const & second = agent.bulk_contains(w.reverse);
    for (std::size_t i = 0; i < union_bins.size(); ++i) union_bins[i] = union_bins[i] || second[i];
    return std::ranges::count(union_bins, true);
}

std::vector<std::size_t> positive_bins(Ibf::membership_agent_type & agent, Window const & w, Options const & o)
{
    auto const & first = agent.bulk_contains(index_key(w, o));
    std::vector<bool> positives(first.begin(), first.end());
    if (o.strand == "both" && w.forward != w.reverse)
    {
        auto const & second = agent.bulk_contains(w.reverse);
        for (std::size_t i = 0; i < positives.size(); ++i) positives[i] = positives[i] || second[i];
    }
    std::vector<std::size_t> bins;
    for (std::size_t i = 0; i < positives.size(); ++i)
        if (positives[i]) bins.push_back(i);
    return bins;
}

Query make_query(std::string id, std::size_t index, std::string const & sequence,
                 Options const & o, Ibf const * ibf, RepeatBed const & query_bed)
{
    Query q;
    q.id = std::move(id);
    q.index = index;
    WindowMaker maker{o, o.repeat == "exclude" && o.mask_target != "reference"};
    BedCursor bed{query_bed};
    bed.reset(q.id);
    std::optional<Ibf::membership_agent_type> agent;
    if (ibf) agent.emplace(ibf->membership_agent());
    q.windows.reserve(sequence.size() >= o.k ? sequence.size() - o.k + 1 : 0);
    for (std::size_t pos = 0; pos < sequence.size(); ++pos)
    {
        auto w = maker.push(sequence[pos], bed.at(pos));
        if (!w) continue;
        ++q.total_windows;
        if (w->base_valid) ++q.base_valid_windows;
        if (w->base_valid && w->repeat) ++q.repeat_windows;
        if (w->valid) ++q.valid_windows;
        if (o.mode == "calibrate")
        {
            q.positive_bins.push_back(w->valid && agent ? positive_bins(*agent, *w, o)
                                                         : std::vector<std::size_t>{});
            q.bin_hits.push_back(q.positive_bins.back().size());
            if (q.bin_hits.back() > o.max_ibf_bin_pairs - q.positive_bin_pairs)
                throw std::runtime_error("IBF-positive bin pairs exceed --max-ibf-bin-pairs for one query; "
                                         "reduce --batch-windows or raise the explicit cap");
            q.positive_bin_pairs += q.bin_hits.back();
        }
        else q.bin_hits.push_back(w->valid && agent ? count_bins(*agent, *w, o) : 0);
        q.windows.push_back(*w);
    }
    return q;
}

ExactResult exact_batch(Options const & o, ReferenceInfo const & info, RepeatBed const & reference_bed,
                        std::vector<Query> const & batch)
{
    std::unordered_set<uint64_t> sought;
    bool collect_hits = o.local_window || !o.positions.empty() || !o.reference_hit_classes.empty() ||
                        !o.exon_blocks.empty() ||
                        !o.loci_output.empty() || !o.expected_loci_bed.empty();
    std::unordered_map<uint64_t, std::vector<std::pair<std::size_t, std::size_t>>> owners;
    ExactResult result;
    if (collect_hits) result.hits.resize(batch.size());
    for (std::size_t qi = 0; qi < batch.size(); ++qi)
        for (std::size_t wi = 0; wi < batch[qi].windows.size(); ++wi)
        {
            auto const & w = batch[qi].windows[wi];
            if (w.valid)
            {
                sought.insert(w.forward);
                if (o.strand != "forward") sought.insert(w.reverse);
                if (collect_hits)
                {
                    owners[w.forward].emplace_back(qi, wi);
                    if (o.strand != "forward" && w.reverse != w.forward)
                        owners[w.reverse].emplace_back(qi, wi);
                }
            }
        }
    if (sought.empty()) return result;
    WindowMaker maker{o, o.repeat == "exclude" && o.mask_target != "query"};
    BedCursor bed{reference_bed};
    std::size_t contig = 0, hit_rows = 0, calibration_pairs = 0;
    scan_fasta(o.reference,
        [&](std::string const & id) { maker.reset(); bed.reset(id); ++contig; },
        [&](std::size_t pos, char c)
        {
            auto w = maker.push(c, bed.at(pos));
            if (!w || !w->valid || !sought.contains(w->forward)) return;
            result.found.insert(w->forward);
            if (o.mode == "calibrate")
            {
                std::size_t bin = info.bin_offsets[contig - 1] + (pos + 1 - o.k) / o.fragment_size;
                auto & bins = result.exact_bins[w->forward];
                if (bins.empty() || bins.back() != bin)
                {
                    bins.push_back(bin);
                    if (++calibration_pairs > o.max_calibration_pairs)
                        throw std::runtime_error("Exact bin pairs exceed --max-calibration-pairs for one batch; "
                                                 "reduce --batch-windows or raise the explicit cap");
                }
            }
            if (!collect_hits) return;
            auto it = owners.find(w->forward);
            if (it == owners.end()) return;
            for (auto const & [qi, wi] : it->second)
            {
                auto const & query_window = batch[qi].windows[wi];
                char strand = query_window.forward == query_window.reverse ? '.' :
                              query_window.forward == w->forward ? '+' : '-';
                result.hits[qi].push_back(Hit{wi, contig - 1, pos + 1 - o.k, strand, w->repeat_class});
                if (++hit_rows > o.max_position_rows)
                    throw std::runtime_error("Exact occurrence count exceeds --max-position-rows for one batch; "
                                             "reduce --batch-windows or raise the explicit cap");
            }
        });
    return result;
}

bool exact_match(Window const & w, Options const & o, std::unordered_set<uint64_t> const & found)
{
    return w.valid && (found.contains(w.forward) ||
           (o.strand != "forward" && found.contains(w.reverse)));
}

void write_fraction(std::ostream & out, std::size_t numerator, std::size_t denominator)
{
    if (denominator) out << std::fixed << std::setprecision(6)
                         << static_cast<double>(numerator) / denominator;
    else out << "NA";
}

struct Totals
{
    std::size_t queries = 0, valid = 0, absent = 0, exact_positive = 0, fp = 0, fn = 0, matched = 0;
    std::size_t exon_eligible = 0, exon_matched = 0;
    std::size_t discovered_locus_rows = 0, expected_valid_windows = 0, expected_locus_matched = 0;
    std::vector<std::size_t> bin_exact_present, bin_fp, bin_fn;
    double query_seconds = 0, exact_seconds = 0, result_seconds = 0;
};

std::vector<std::size_t> exact_bins_for_window(Window const & w, Options const & o,
                                                ExactResult const & result)
{
    std::vector<std::size_t> bins;
    auto append = [&](uint64_t key)
    {
        auto it = result.exact_bins.find(key);
        if (it != result.exact_bins.end()) bins.insert(bins.end(), it->second.begin(), it->second.end());
    };
    append(w.forward);
    if (o.strand != "forward" && w.reverse != w.forward) append(w.reverse);
    std::sort(bins.begin(), bins.end());
    bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
    return bins;
}

struct LocalStats
{
    std::size_t best = 0, regions = 0, contig = 0, start = 0, end = 0, multimapped = 0;
    std::size_t sliding_best = 0, sliding_contig = 0, sliding_start = 0, sliding_end = 0;
};

LocalStats local_stats(std::vector<Hit> const & hits, std::size_t window_count,
                       ReferenceInfo const & info, std::size_t region_length)
{
    LocalStats stats;
    std::vector<unsigned char> occurrence_count(window_count, 0);
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> region_windows;
    region_windows.reserve(hits.size());
    for (auto const & hit : hits)
    {
        if (occurrence_count[hit.window] == 1) ++stats.multimapped;
        if (occurrence_count[hit.window] < 2) ++occurrence_count[hit.window];
        region_windows.emplace_back(hit.contig, hit.start / region_length, hit.window);
    }
    std::sort(region_windows.begin(), region_windows.end());
    region_windows.erase(std::unique(region_windows.begin(), region_windows.end()), region_windows.end());
    std::size_t i = 0;
    while (i < region_windows.size())
    {
        auto [contig, region, ignored_window] = region_windows[i];
        (void) ignored_window;
        std::size_t j = i + 1;
        while (j < region_windows.size() && std::get<0>(region_windows[j]) == contig &&
               std::get<1>(region_windows[j]) == region) ++j;
        ++stats.regions;
        if (j - i > stats.best)
        {
            stats.best = j - i;
            stats.contig = contig;
            stats.start = region * region_length;
            stats.end = stats.start + std::min(region_length, info.lengths[contig] - stats.start);
        }
        i = j;
    }
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> events;
    events.reserve(hits.size());
    for (auto const & hit : hits) events.emplace_back(hit.contig, hit.start, hit.window);
    std::sort(events.begin(), events.end());
    events.erase(std::unique(events.begin(), events.end()), events.end());
    std::unordered_map<std::size_t, std::size_t> window_frequency;
    std::size_t right = 0;
    for (std::size_t left = 0; left < events.size(); ++left)
    {
        auto const [contig, start, ignored] = events[left];
        (void) ignored;
        while (right < events.size() && std::get<0>(events[right]) == contig &&
               std::get<1>(events[right]) - start < region_length)
            ++window_frequency[std::get<2>(events[right++])];
        if (window_frequency.size() > stats.sliding_best)
        {
            stats.sliding_best = window_frequency.size();
            stats.sliding_contig = contig;
            stats.sliding_start = start;
            stats.sliding_end = start + std::min(region_length, info.lengths[contig] - start);
        }
        auto it = window_frequency.find(std::get<2>(events[left]));
        if (--it->second == 0) window_frequency.erase(it);
    }
    return stats;
}

struct LocusCandidate
{
    std::size_t contig = 0, start = 0, end = 0, matched_windows = 0, occurrence_count = 0;
};

std::vector<LocusCandidate> discover_loci(std::vector<Hit> const & hits, Options const & o)
{
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> events;
    events.reserve(hits.size());
    for (auto const & hit : hits) events.emplace_back(hit.contig, hit.start, hit.window);
    std::sort(events.begin(), events.end());
    events.erase(std::unique(events.begin(), events.end()), events.end());
    std::vector<LocusCandidate> candidates;
    std::unordered_set<std::size_t> windows;
    LocusCandidate current;
    std::size_t previous_start = 0;
    bool active = false;
    auto finish = [&]()
    {
        if (!active) return;
        current.matched_windows = windows.size();
        candidates.push_back(current);
        windows.clear();
    };
    std::size_t effective_gap = std::max(o.locus_gap, o.k - 1);
    for (auto const & [contig, start, window] : events)
    {
        bool connected = active && contig == current.contig && start - previous_start <= effective_gap;
        if (connected && start + o.k - current.start > o.max_locus_span)
        {
            if (start < current.end)
                throw std::runtime_error("Overlapping exact-hit chain exceeds --max-locus-span; "
                                         "raise the cap or reduce --locus-gap");
            connected = false;
        }
        if (!connected)
        {
            finish();
            current = {contig, start, start + o.k, 0, 0};
            active = true;
        }
        current.end = std::max(current.end, start + o.k);
        ++current.occurrence_count;
        windows.insert(window);
        previous_start = start;
    }
    finish();
    std::sort(candidates.begin(), candidates.end(), [](auto const & a, auto const & b)
    {
        if (a.matched_windows != b.matched_windows) return a.matched_windows > b.matched_windows;
        return std::tie(a.contig, a.start, a.end) < std::tie(b.contig, b.start, b.end);
    });
    return candidates;
}

struct ExpectedStats
{
    std::size_t matched_windows = 0, hit_rows = 0;
    std::vector<std::size_t> region_matched, region_hits;
};

ExpectedStats expected_locus_stats(std::vector<Hit> const & hits, std::size_t query_windows,
                                   std::vector<ExpectedLocus> const & regions, std::size_t k)
{
    ExpectedStats stats;
    stats.region_matched.resize(regions.size());
    stats.region_hits.resize(regions.size());
    std::vector<unsigned char> matched(query_windows, 0);
    std::vector<std::unordered_set<std::size_t>> region_windows(regions.size());
    for (auto const & hit : hits)
    {
        auto it = std::upper_bound(regions.begin(), regions.end(), hit,
            [](Hit const & h, ExpectedLocus const & region)
            { return std::tie(h.contig, h.start) < std::tie(region.contig, region.start); });
        if (it == regions.begin()) continue;
        --it;
        if (it->contig != hit.contig || hit.start < it->start || hit.start + k > it->end) continue;
        std::size_t region = static_cast<std::size_t>(it - regions.begin());
        ++stats.region_hits[region];
        ++stats.hit_rows;
        if (!matched[hit.window]) { matched[hit.window] = 1; ++stats.matched_windows; }
        region_windows[region].insert(hit.window);
    }
    for (std::size_t i = 0; i < regions.size(); ++i) stats.region_matched[i] = region_windows[i].size();
    return stats;
}

void write_batch(std::ostream & out, Options const & o, std::string const & ref_id,
                 std::vector<Query> const & batch, ExactResult const & exact_result,
                 ReferenceInfo const & info, PositionWriter & position_writer, Totals & totals,
                 std::ostream * repeat_class_out, RepeatBed const & query_bed, ExonBlocks const & exon_blocks,
                 std::ostream * reference_hit_class_out, RepeatBed const & reference_bed,
                 std::ostream * loci_out, std::ostream * expected_out, ExpectedLoci const & expected_loci)
{
    for (std::size_t qi = 0; qi < batch.size(); ++qi)
    {
        auto const & q = batch[qi];
        std::size_t ibf_match = 0, exact_count = 0, fp = 0, fn = 0, bin_hits = 0;
        std::size_t repeat_match = 0, nonrepeat_valid = 0, nonrepeat_match = 0;
        std::map<uint32_t, std::pair<std::size_t, std::size_t>> class_counts;
        std::unordered_set<uint64_t> distinct, distinct_match;
        for (std::size_t i = 0; i < q.windows.size(); ++i)
        {
            auto const & w = q.windows[i];
            if (!w.valid) continue;
            if (o.mode == "calibrate")
            {
                auto exact_bins = exact_bins_for_window(w, o, exact_result);
                auto const & positive = q.positive_bins[i];
                for (std::size_t bin : exact_bins)
                {
                    ++totals.bin_exact_present[bin];
                    if (!std::binary_search(positive.begin(), positive.end(), bin)) ++totals.bin_fn[bin];
                }
                for (std::size_t bin : positive)
                    if (!std::binary_search(exact_bins.begin(), exact_bins.end(), bin)) ++totals.bin_fp[bin];
            }
            uint64_t key = std::min(w.forward, w.reverse);
            distinct.insert(key);
            bool ibf = q.bin_hits[i] > 0;
            bool exact = o.mode != "ibf" && exact_match(w, o, exact_result.found);
            if (ibf) ++ibf_match;
            if (exact) ++exact_count;
            if (o.mode == "ibf" ? ibf : exact) distinct_match.insert(key);
            if (o.mode != "ibf" && ibf && !exact) ++fp;
            if (o.mode != "ibf" && !ibf && exact && o.mode != "exact") ++fn;
            bin_hits += q.bin_hits[i];
            bool selected_match = o.mode == "ibf" ? ibf : exact;
            if (repeat_class_out && w.repeat_class)
            {
                auto & counts = class_counts[w.repeat_class];
                ++counts.first;
                if (selected_match) ++counts.second;
            }
            if (!w.repeat) { ++nonrepeat_valid; if (selected_match) ++nonrepeat_match; }
            else if (selected_match) ++repeat_match;
        }
        std::size_t matched = o.mode == "ibf" ? ibf_match : exact_count;
        std::string status = q.valid_windows == 0 ? "NO_VALID_KMERS" :
                             matched == 0 ? "NO_MATCH" :
                             static_cast<double>(matched) / q.valid_windows >= o.strict ? "PASS_STRICT" :
                             static_cast<double>(matched) / q.valid_windows >= o.basic ? "PASS_BASIC" : "LOW_SUPPORT";
        out << "2\t" << q.index << '\t' << q.id << '\t' << ref_id << '\t' << o.k << '\t'
            << q.total_windows << '\t' << q.valid_windows << '\t'
            << q.total_windows - q.valid_windows << '\t' << matched << '\t'
            << q.valid_windows - matched << '\t' << distinct.size() << '\t'
            << distinct_match.size() << '\t' << bin_hits << '\t';
        write_fraction(out, matched, q.valid_windows);
        out << '\t' << status << '\t' << o.strand << '\t' << o.repeat << '\t' << o.ambiguity
            << '\t' << o.mode << '\t' << o.fragment_size << '\t' << o.hash_functions << '\t'
            << std::fixed << std::setprecision(o.fpr < 0.000001 ? 12 : 6) << o.fpr << '\t'
            << std::setprecision(6) << o.basic << '\t' << o.strict
            << '\t' << o.threshold_source << '\t'
            << (o.training_control_set.empty() ? "NA" : o.training_control_set) << '\t'
            << (o.evaluation_control_set.empty() ? "NA" : o.evaluation_control_set) << '\t';
        if (o.mode == "exact") out << "NA"; else write_fraction(out, ibf_match, q.valid_windows);
        out << '\t';
        if (o.mode == "ibf") out << "NA"; else write_fraction(out, exact_count, q.valid_windows);
        out << '\t';
        if (o.mode == "ibf" || o.mode == "exact") out << "NA\tNA\tNA";
        else
        {
            out << fp << '\t' << fn << '\t';
            if (q.valid_windows) out << std::fixed << std::setprecision(6)
                                     << (static_cast<double>(ibf_match) - exact_count) / q.valid_windows;
            else out << "NA";
        }
        out << '\t'; write_fraction(out, q.repeat_windows, q.base_valid_windows);
        out << '\t' << nonrepeat_valid << '\t' << nonrepeat_match << '\t';
        write_fraction(out, nonrepeat_match, nonrepeat_valid);
        out << '\t' << repeat_match << '\t';
        std::vector<Hit> const * hits = exact_result.hits.empty() ? nullptr : &exact_result.hits[qi];
        if (reference_hit_class_out)
        {
            std::map<uint32_t, std::pair<std::unordered_set<std::size_t>, std::size_t>> class_hits;
            for (Hit const & hit : *hits)
            {
                auto & counts = class_hits[hit.reference_repeat_class];
                counts.first.insert(hit.window);
                ++counts.second;
            }
            for (auto const & [class_id, counts] : class_hits)
            {
                std::string label = class_id == mixed_repeat_class ? "MULTIPLE_CLASSES" :
                                    class_id < reference_bed.classes.size() ? reference_bed.classes[class_id] :
                                    "UNKNOWN";
                *reference_hit_class_out << "2\t" << q.index << '\t' << q.id << '\t' << ref_id << '\t'
                                         << label << '\t' << counts.first.size() << '\t' << counts.second
                                         << '\t' << exact_count << '\t';
                write_fraction(*reference_hit_class_out, counts.first.size(), exact_count);
                *reference_hit_class_out << '\n';
            }
        }
        if (o.local_window)
        {
            LocalStats local = local_stats(*hits, q.windows.size(), info, o.local_window);
            write_fraction(out, exact_count, q.valid_windows);
            out << '\t'; write_fraction(out, local.best, q.valid_windows);
            out << '\t' << local.best << '\t' << local.regions << '\t';
            if (local.regions)
                out << info.ids[local.contig] << '\t' << local.start << '\t' << local.end;
            else out << "NA\tNA\tNA";
            out << '\t'; write_fraction(out, local.best, exact_count);
            out << '\t';
            if (exact_count) out << std::fixed << std::setprecision(6)
                                 << 1.0 - static_cast<double>(local.best) / exact_count;
            else out << "NA";
            out << '\t' << local.multimapped;
            out << '\t'; write_fraction(out, local.sliding_best, q.valid_windows);
            out << '\t' << local.sliding_best << '\t';
            if (local.sliding_best)
                out << info.ids[local.sliding_contig] << '\t' << local.sliding_start << '\t'
                    << local.sliding_end;
            else out << "NA\tNA\tNA";
            out << '\t' << local.sliding_best - local.best;
        }
        else out << "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
        auto block_it = exon_blocks.find(q.id);
        if (!o.exon_blocks.empty() && block_it != exon_blocks.end())
        {
            std::size_t eligible = 0, exon_matched = 0, block_index = 0;
            std::unordered_multimap<std::size_t, Hit const *> hits_by_window;
            if (hits) for (auto const & hit : *hits) hits_by_window.emplace(hit.window, &hit);
            auto const & blocks = block_it->second;
            for (std::size_t wi = 0; wi < q.windows.size(); ++wi)
            {
                if (!q.windows[wi].valid) continue;
                while (block_index < blocks.size() && blocks[block_index].query_end <= wi) ++block_index;
                if (block_index == blocks.size()) break;
                auto const & block = blocks[block_index];
                if (wi < block.query_start || wi + o.k > block.query_end) continue;
                ++eligible;
                std::size_t expected = block.strand == '+' ? block.ref_start + wi - block.query_start :
                                       block.ref_end - (wi - block.query_start) - o.k;
                auto [first_hit, last_hit] = hits_by_window.equal_range(wi);
                for (auto it = first_hit; it != last_hit; ++it)
                {
                    Hit const & hit = *it->second;
                    if (hit.contig == block.contig && hit.start == expected &&
                        (hit.strand == block.strand || hit.strand == '.'))
                    {
                        ++exon_matched;
                        break;
                    }
                }
            }
            out << '\t' << eligible << '\t' << exon_matched << '\t';
            write_fraction(out, exon_matched, eligible);
            totals.exon_eligible += eligible;
            totals.exon_matched += exon_matched;
        }
        else out << "\tNA\tNA\tNA";
        if (loci_out)
        {
            auto candidates = discover_loci(*hits, o);
            if (candidates.size() > o.max_locus_rows - totals.discovered_locus_rows)
                throw std::runtime_error("Discovered loci exceed --max-locus-rows");
            totals.discovered_locus_rows += candidates.size();
            out << '\t' << candidates.size() << '\t';
            if (candidates.empty()) out << "0\tNA\tNA\tNA\tNA";
            else
            {
                auto const & best = candidates.front();
                out << best.matched_windows << '\t';
                write_fraction(out, best.matched_windows, q.valid_windows);
                out << '\t' << info.ids[best.contig] << '\t' << best.start << '\t' << best.end;
            }
            for (std::size_t ci = 0; ci < candidates.size(); ++ci)
            {
                auto const & locus = candidates[ci];
                *loci_out << "2\t" << q.index << '\t' << q.id << '\t' << ref_id << '\t' << ci + 1
                          << '\t' << info.ids[locus.contig] << '\t' << locus.start << '\t' << locus.end
                          << '\t' << locus.matched_windows << '\t' << locus.occurrence_count << '\t'
                          << q.valid_windows << '\t';
                write_fraction(*loci_out, locus.matched_windows, q.valid_windows);
                *loci_out << '\t'; write_fraction(*loci_out, locus.matched_windows, exact_count);
                *loci_out << '\n';
            }
        }
        else out << "\tNA\tNA\tNA\tNA\tNA\tNA";
        auto expected_it = expected_loci.find(q.id);
        if (expected_out && expected_it != expected_loci.end())
        {
            auto const & regions = expected_it->second;
            ExpectedStats stats = expected_locus_stats(*hits, q.windows.size(), regions, o.k);
            out << '\t' << regions.size() << '\t' << stats.matched_windows << '\t';
            write_fraction(out, stats.matched_windows, q.valid_windows);
            out << '\t' << exact_count - stats.matched_windows << '\t' << stats.hit_rows;
            totals.expected_valid_windows += q.valid_windows;
            totals.expected_locus_matched += stats.matched_windows;
            for (std::size_t ri = 0; ri < regions.size(); ++ri)
            {
                auto const & region = regions[ri];
                *expected_out << "2\t" << q.index << '\t' << q.id << '\t' << ref_id << '\t'
                              << info.ids[region.contig] << '\t' << region.start << '\t' << region.end
                              << '\t' << region.label << '\t' << q.valid_windows << '\t'
                              << stats.region_matched[ri] << '\t';
                write_fraction(*expected_out, stats.region_matched[ri], q.valid_windows);
                *expected_out << '\t' << stats.region_hits[ri] << '\n';
            }
        }
        else out << "\tNA\tNA\tNA\tNA\tNA";
        out << '\n';
        if (repeat_class_out)
            for (auto const & [class_id, counts] : class_counts)
            {
                std::string label = class_id == mixed_repeat_class ? "MULTIPLE_CLASSES" :
                                    class_id < query_bed.classes.size() ? query_bed.classes[class_id] : "UNKNOWN";
                *repeat_class_out << "2\t" << q.index << '\t' << q.id << '\t' << label << '\t'
                                  << counts.first << '\t' << counts.second << '\t';
                write_fraction(*repeat_class_out, counts.second, counts.first);
                *repeat_class_out << '\n';
            }
        if (!o.positions.empty())
            for (auto const & hit : *hits)
                position_writer.write(std::to_string(q.index) + '\t' + q.id + '\t' + ref_id + '\t' +
                                      info.ids[hit.contig] + '\t' + std::to_string(hit.start) + '\t' +
                                      std::to_string(hit.start + o.k) + '\t' + hit.strand + '\t' +
                                      std::to_string(hit.window + 1) + '\n');
        if (fn) std::cerr << "Diagnostic: " << q.id << " has " << fn
                          << " exact-positive/IBF-negative windows; check preprocessing, strand policy and index version.\n";
        ++totals.queries;
        totals.valid += q.valid_windows;
        totals.matched += matched;
        totals.fp += fp;
        totals.fn += fn;
        if (o.mode != "ibf")
        {
            totals.absent += q.valid_windows - exact_count;
            totals.exact_positive += exact_count;
        }
    }
}

std::size_t peak_rss_bytes()
{
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
#if defined(__APPLE__)
    return static_cast<std::size_t>(usage.ru_maxrss);
#else
    return static_cast<std::size_t>(usage.ru_maxrss) * 1024;
#endif
}

double elapsed(Clock::time_point start) { return std::chrono::duration<double>(Clock::now() - start).count(); }

void write_fragments(fs::path const & path, Options const & o, ReferenceInfo const & info, bool indexed)
{
    if (path.empty()) return;
    fs::create_directories(path.parent_path().empty() ? fs::path{"."} : path.parent_path());
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write " + path.string());
    out << "schema_version\tbin_index\tcontig\tstart\tend\tsequence_end\tinserted_reference_windows"
           "\tdistinct_reference_kmers\n";
    for (std::size_t ci = 0; ci < info.ids.size(); ++ci)
        for (std::size_t bin = info.bin_offsets[ci];
             bin < info.bin_offsets[ci] + (info.lengths[ci] + o.fragment_size - 1) / o.fragment_size; ++bin)
        {
            std::size_t start = (bin - info.bin_offsets[ci]) * o.fragment_size;
            std::size_t end = std::min(start + o.fragment_size, info.lengths[ci]);
            std::size_t sequence_end = std::min(end + o.k - 1, info.lengths[ci]);
            out << "2\t" << bin << '\t' << info.ids[ci] << '\t' << start << '\t' << end << '\t'
                << sequence_end << '\t';
            if (indexed) out << info.bin_insertions[bin] << '\t' << info.bin_distinct[bin];
            else out << "NA\tNA";
            out << '\n';
        }
    if (!out) throw std::runtime_error("Failed writing fragment coordinates");
}

double observed_ibf_union_estimate(Ibf const & ibf, ReferenceInfo const & info, Options const & o)
{
    std::size_t technical_bins = ibf.raw_data().size() / info.bin_bits;
    double estimate = 0;
    for (std::size_t bin = 0; bin < info.bin_count; ++bin)
    {
        std::size_t occupied = 0;
        for (std::size_t bit = 0; bit < info.bin_bits; ++bit)
            occupied += ibf.raw_data()[bit * technical_bins + bin];
        estimate += std::pow(static_cast<double>(occupied) / info.bin_bits,
                             static_cast<double>(o.hash_functions));
    }
    return std::min(1.0, estimate);
}

int run(Options o)
{
    auto run_start = Clock::now();
    auto reference_start = Clock::now();
    std::string reference_bed_sha = sha256_file(o.reference_repeat_bed);
    std::string query_bed_sha = sha256_file(o.query_repeat_bed);
    std::string exon_blocks_sha = sha256_file(o.exon_blocks);
    std::string expected_loci_sha = sha256_file(o.expected_loci_bed);
    std::string threshold_table_sha = sha256_file(o.threshold_table);
    std::string config_sha = sha256_file(o.config_path);
    RepeatBed reference_bed = load_repeat_bed(o.reference_repeat_bed, o.max_repeat_intervals);
    RepeatBed query_bed = load_repeat_bed(o.query_repeat_bed, o.max_repeat_intervals);
    ReferenceInfo info = inspect_reference(o, reference_bed);
    if (o.mode == "ibf" || o.genome_fpr_given)
        o.fpr = std::min(o.fpr, o.genome_fpr / std::max<std::size_t>(1, info.bin_count));
    double nominal_any_bin_fpr = -std::expm1(static_cast<double>(info.bin_count) * std::log1p(-o.fpr));
    if (o.mode != "exact" && nominal_any_bin_fpr > 0.1)
        std::cerr << "Warning: nominal any-bin false-positive probability is " << nominal_any_bin_fpr
                  << " across " << info.bin_count << " bins. IBF-only recovery can be unreliable; "
                     "use exact-verified results for classification.\n";
    for (auto const & [id, intervals] : reference_bed.intervals)
    {
        auto it = std::find(info.ids.begin(), info.ids.end(), id);
        if (it == info.ids.end()) throw std::runtime_error("Repeat BED contig missing from reference: " + id);
        std::size_t length = info.lengths[static_cast<std::size_t>(it - info.ids.begin())];
        if (!intervals.empty() && intervals.back().start >= length)
            throw std::runtime_error("Repeat BED interval outside reference contig: " + id);
        for (auto const & interval : intervals)
            if (interval.end > length)
                throw std::runtime_error("Repeat BED interval outside reference contig: " + id);
    }
    double inspect_seconds = elapsed(reference_start);
    std::size_t inspect_peak = peak_rss_bytes();
    auto index_start = Clock::now();
    std::unique_ptr<Ibf> ibf;
    double index_seconds = 0, index_load_seconds = 0, index_save_seconds = 0;
    std::size_t disk_index_bytes = 0;
    std::string index_sha, index_provenance = "none";
    fs::path index_path = !o.load_index.empty() ? o.load_index : o.save_index;
    if (!o.load_index.empty())
    {
        IndexMetadata expected = make_index_metadata(o, info, reference_bed_sha);
        ibf = load_index_file(o.load_index, expected, info);
        disk_index_bytes = fs::file_size(o.load_index);
        index_sha = sha256_file(o.load_index);
        index_load_seconds = elapsed(index_start);
        index_provenance = "loaded";
    }
    else if (o.mode != "exact")
    {
        ibf = build_index(o, info, reference_bed);
        index_seconds = elapsed(index_start);
        index_provenance = "built";
        if (!o.save_index.empty())
        {
            auto save_start = Clock::now();
            IndexMetadata meta = make_index_metadata(o, info, reference_bed_sha);
            save_index_file(o.save_index, meta, *ibf);
            disk_index_bytes = fs::file_size(o.save_index);
            index_sha = sha256_file(o.save_index);
            write_index_metadata_json(o.save_index, meta, index_sha, disk_index_bytes);
            index_save_seconds = elapsed(save_start);
        }
    }
    std::size_t index_peak = peak_rss_bytes();
    double occupancy_estimate = -1, occupancy_guard_seconds = 0;
    std::size_t occupancy_guard_peak = 0;
    if (o.mode == "ibf")
    {
        auto guard_start = Clock::now();
        occupancy_estimate = observed_ibf_union_estimate(*ibf, info, o);
        occupancy_guard_seconds = elapsed(guard_start);
        occupancy_guard_peak = peak_rss_bytes();
        if (occupancy_estimate > o.genome_fpr * (1.0 + 1e-9))
            throw std::runtime_error("IBF occupancy-based genome-wide false-positive estimate exceeds "
                                     "--genome-fpr; use exact mode or rebuild with a smaller --fpr");
    }
    auto fragment_start = Clock::now();
    write_fragments(o.fragments, o, info, ibf != nullptr);
    double fragment_seconds = elapsed(fragment_start);
    std::size_t fragment_peak = peak_rss_bytes();

    fs::create_directories(o.output.parent_path().empty() ? fs::path{"."} : o.output.parent_path());
    std::ofstream out(o.output);
    if (!out) throw std::runtime_error("Cannot write " + o.output.string());
    out << "schema_version\tquery_index\tquery_id\treference_id\tk\ttotal_windows\tvalid_windows\tinvalid_windows"
           "\tmatched_windows\tunmatched_windows\tdistinct_query_kmers\tdistinct_matched_kmers\ttotal_bin_hits"
           "\trecovery_fraction\tstatus\tstrand_mode\trepeat_policy\tambiguity_policy\tmode"
           "\tfragment_size\thash_functions\tconfigured_per_bin_fpr\tbasic_threshold\tstrict_threshold"
           "\tthreshold_source\ttraining_control_set\tevaluation_control_set"
           "\tibf_recovery"
           "\texact_recovery\tfalse_positive_windows\tfalse_negative_windows\trecovery_difference"
           "\trepeat_fraction\tvalid_nonrepeat_windows\tmatched_nonrepeat_windows\tnonrepeat_recovery"
           "\trepeat_matched_windows\tglobal_recovery\tbest_region_recovery"
           "\tbest_region_matched_windows\tnumber_of_supported_regions\tcandidate_contig"
           "\tcandidate_start\tcandidate_end\tlocal_concentration\toff_region_fraction"
           "\tmultimapped_query_windows\tsliding_best_region_recovery\tsliding_best_region_matched_windows"
           "\tsliding_candidate_contig\tsliding_candidate_start\tsliding_candidate_end"
           "\tboundary_gain_windows\texon_eligible_windows\texon_matched_windows"
           "\texon_aware_recovery\tnumber_of_discovered_loci\tbest_locus_matched_windows"
           "\tbest_locus_recovery\tbest_locus_contig\tbest_locus_start\tbest_locus_end"
           "\texpected_locus_count\texpected_locus_matched_windows\texpected_locus_recovery"
           "\toff_expected_locus_windows\texpected_locus_hit_rows\n";
    PositionWriter position_writer{o.positions};
    std::ofstream repeat_class_file;
    if (!o.repeat_classes.empty())
    {
        fs::create_directories(o.repeat_classes.parent_path().empty() ? fs::path{"."} :
                               o.repeat_classes.parent_path());
        repeat_class_file.open(o.repeat_classes);
        if (!repeat_class_file) throw std::runtime_error("Cannot write " + o.repeat_classes.string());
        repeat_class_file << "schema_version\tquery_index\tquery_id\trepeat_class\tvalid_windows"
                             "\tmatched_windows\trecovery_fraction\n";
    }
    std::ofstream reference_hit_class_file;
    if (!o.reference_hit_classes.empty())
    {
        fs::create_directories(o.reference_hit_classes.parent_path().empty() ? fs::path{"."} :
                               o.reference_hit_classes.parent_path());
        reference_hit_class_file.open(o.reference_hit_classes);
        if (!reference_hit_class_file)
            throw std::runtime_error("Cannot write " + o.reference_hit_classes.string());
        reference_hit_class_file << "schema_version\tquery_index\tquery_id\treference_id"
                                    "\treference_repeat_class\tmatched_query_windows\texact_occurrences"
                                    "\texact_matched_query_windows\tfraction_of_exact_matched_windows\n";
    }
    std::ofstream loci_file;
    if (!o.loci_output.empty())
    {
        fs::create_directories(o.loci_output.parent_path().empty() ? fs::path{"."} : o.loci_output.parent_path());
        loci_file.open(o.loci_output);
        if (!loci_file) throw std::runtime_error("Cannot write " + o.loci_output.string());
        loci_file << "schema_version\tquery_index\tquery_id\treference_id\trank\tcontig\tstart\tend"
                     "\tmatched_windows\texact_occurrences\tvalid_query_windows\tlocus_recovery"
                     "\tlocus_concentration\n";
    }
    std::ofstream expected_file;
    if (!o.expected_loci_output.empty())
    {
        fs::create_directories(o.expected_loci_output.parent_path().empty() ? fs::path{"."} :
                               o.expected_loci_output.parent_path());
        expected_file.open(o.expected_loci_output);
        if (!expected_file) throw std::runtime_error("Cannot write " + o.expected_loci_output.string());
        expected_file << "schema_version\tquery_index\tquery_id\treference_id\tcontig\tstart\tend"
                         "\tlabel\tvalid_query_windows\tmatched_windows\tlocus_recovery\texact_occurrences\n";
    }
    std::string ref_id = reference_id_from_path(o.reference);
    ExonBlocks exon_blocks = load_exon_blocks(o, info);
    ExpectedLoci expected_loci = load_expected_loci(o, info);

    Totals totals;
    if (o.mode == "calibrate")
    {
        totals.bin_exact_present.resize(info.bin_count);
        totals.bin_fp.resize(info.bin_count);
        totals.bin_fn.resize(info.bin_count);
    }
    std::vector<Query> batch;
    std::size_t batch_size = 0, query_index = 0;
    std::size_t batch_positive_bin_pairs = 0;
    std::string query_id, sequence;
    std::unordered_set<std::string> seen_query_ids;
    uint64_t query_sequence_hash = 14695981039346656037ULL;
    fraglncrna::Sha256 query_sha;
    query_sha.update("FragIncRNA-FASTA-semantic-v1");
    std::size_t query_sequence_length = 0;
    bool first = true;
    auto flush = [&]()
    {
        if (batch.empty()) return;
        auto t = Clock::now();
        ExactResult exact_result;
        if (o.mode != "ibf") exact_result = exact_batch(o, info, reference_bed, batch);
        totals.exact_seconds += elapsed(t);
        t = Clock::now();
        write_batch(out, o, ref_id, batch, exact_result, info, position_writer, totals,
                    repeat_class_file.is_open() ? &repeat_class_file : nullptr, query_bed, exon_blocks,
                    reference_hit_class_file.is_open() ? &reference_hit_class_file : nullptr, reference_bed,
                    loci_file.is_open() ? &loci_file : nullptr,
                    expected_file.is_open() ? &expected_file : nullptr, expected_loci);
        totals.result_seconds += elapsed(t);
        batch.clear(); batch_size = 0; batch_positive_bin_pairs = 0;
    };
    auto finish_query = [&]()
    {
        if (first) return;
        auto query_intervals = query_bed.intervals.find(query_id);
        if (query_intervals != query_bed.intervals.end())
            for (auto const & interval : query_intervals->second)
                if (interval.end > sequence.size())
                    throw std::runtime_error("Repeat BED interval outside query: " + query_id);
        auto exon_it = exon_blocks.find(query_id);
        if (exon_it != exon_blocks.end())
            for (auto const & block : exon_it->second)
                if (block.query_end > sequence.size())
                    throw std::runtime_error("Exon block outside query: " + query_id);
        auto t = Clock::now();
        Query q = make_query(query_id, query_index++, sequence, o, ibf.get(), query_bed);
        totals.query_seconds += elapsed(t);
        if (o.mode == "calibrate")
        {
            if (q.positive_bin_pairs > o.max_ibf_bin_pairs - batch_positive_bin_pairs)
                flush();
            batch_positive_bin_pairs += q.positive_bin_pairs;
        }
        batch_size += q.windows.size();
        batch.push_back(std::move(q));
        sequence.clear();
        if (batch_size >= o.batch_windows) flush();
    };
    scan_fasta(o.queries,
        [&](std::string const & id)
        {
            if (!seen_query_ids.insert(id).second)
                throw std::runtime_error("Duplicate query FASTA identifier: " + id);
            if (!first) sha256_record_end(query_sha, query_sequence_length);
            finish_query(); query_id = id; first = false;
            query_sequence_length = 0;
            query_sequence_hash = fnv_update(query_sequence_hash, 255);
            sha256_record_start(query_sha, id);
            for (unsigned char c : id) query_sequence_hash = fnv_update(query_sequence_hash, c);
        },
        [&](std::size_t, char c)
        {
            if (sequence.size() >= o.max_query_bases)
                throw std::runtime_error("Query record exceeds --max-query-bases: " + query_id);
            sequence += c;
            ++query_sequence_length;
            query_sequence_hash = fnv_update(query_sequence_hash, static_cast<unsigned char>(c));
            query_sha.update(static_cast<unsigned char>(c));
        });
    if (!first) sha256_record_end(query_sha, query_sequence_length);
    finish_query();
    for (auto const & [id, ignored] : query_bed.intervals)
        if (!seen_query_ids.contains(id)) throw std::runtime_error("Repeat BED query missing from FASTA: " + id);
    for (auto const & [id, ignored] : exon_blocks)
        if (!seen_query_ids.contains(id)) throw std::runtime_error("Exon-block query missing from FASTA: " + id);
    for (auto const & [id, ignored] : expected_loci)
        if (!seen_query_ids.contains(id)) throw std::runtime_error("Expected-locus query missing from FASTA: " + id);
    flush();
    position_writer.close();
    if (repeat_class_file.is_open())
    {
        repeat_class_file.close();
        if (!repeat_class_file) throw std::runtime_error("Failed writing repeat classes");
    }
    if (reference_hit_class_file.is_open())
    {
        reference_hit_class_file.close();
        if (!reference_hit_class_file) throw std::runtime_error("Failed writing reference hit classes");
    }
    if (loci_file.is_open())
    {
        loci_file.close();
        if (!loci_file) throw std::runtime_error("Failed writing discovered loci");
    }
    if (expected_file.is_open())
    {
        expected_file.close();
        if (!expected_file) throw std::runtime_error("Failed writing expected loci");
    }
    if (!out) throw std::runtime_error("Failed writing result TSV");
    out.close();
    std::size_t absent_bin_pairs = 0, fp_bin_pairs = 0, fn_bin_pairs = 0;
    std::size_t occupied_bits_total = 0;
    double occupancy_seconds = 0;
    std::size_t occupancy_peak = 0;
    if (o.mode == "calibrate")
    {
        auto occupancy_start = Clock::now();
        std::vector<std::size_t> occupied(info.bin_count, 0);
        std::size_t technical_bins = ibf->raw_data().size() / info.bin_bits;
        for (std::size_t bit = 0; bit < info.bin_bits; ++bit)
            for (std::size_t bin = 0; bin < info.bin_count; ++bin)
                occupied[bin] += ibf->raw_data()[bit * technical_bins + bin];
        occupied_bits_total = std::accumulate(occupied.begin(), occupied.end(), std::size_t{0});
        occupancy_seconds = elapsed(occupancy_start);
        occupancy_peak = peak_rss_bytes();
        fs::create_directories(o.calibration_bins.parent_path().empty() ? fs::path{"."} :
                               o.calibration_bins.parent_path());
        std::ofstream bins_out(o.calibration_bins);
        if (!bins_out) throw std::runtime_error("Cannot write " + o.calibration_bins.string());
        bins_out << "schema_version\tbin_index\tcontig\tstart\tend\tinserted_reference_windows"
                    "\tdistinct_reference_kmers\tbin_bits\toccupied_bits\tbit_occupancy_fraction"
                    "\tinsertions_per_bit\ttested_absent_query_windows\tfalse_positive_bin_pairs"
                    "\tfalse_negative_bin_pairs\tempirical_per_bin_fpr\n";
        for (std::size_t ci = 0; ci < info.ids.size(); ++ci)
            for (std::size_t bin = info.bin_offsets[ci];
                 bin < info.bin_offsets[ci] + (info.lengths[ci] + o.fragment_size - 1) / o.fragment_size; ++bin)
            {
                std::size_t start = (bin - info.bin_offsets[ci]) * o.fragment_size;
                std::size_t absent = totals.valid - totals.bin_exact_present[bin];
                absent_bin_pairs += absent;
                fp_bin_pairs += totals.bin_fp[bin];
                fn_bin_pairs += totals.bin_fn[bin];
                bins_out << "2\t" << bin << '\t' << info.ids[ci] << '\t' << start << '\t'
                         << std::min(start + o.fragment_size, info.lengths[ci]) << '\t'
                         << info.bin_insertions[bin] << '\t' << info.bin_distinct[bin] << '\t'
                         << info.bin_bits << '\t' << occupied[bin] << '\t'
                         << std::fixed << std::setprecision(6)
                         << static_cast<double>(occupied[bin]) / info.bin_bits << '\t'
                         << static_cast<double>(info.bin_insertions[bin]) / info.bin_bits << '\t'
                         << absent << '\t' << totals.bin_fp[bin] << '\t' << totals.bin_fn[bin] << '\t';
                write_fraction(bins_out, totals.bin_fp[bin], absent);
                bins_out << '\n';
            }
        if (!bins_out) throw std::runtime_error("Failed writing calibration bins");
    }
    fs::create_directories(o.metrics.parent_path().empty() ? fs::path{"."} : o.metrics.parent_path());
    std::ofstream metrics(o.metrics);
    if (!metrics) throw std::runtime_error("Cannot write metrics JSON");
    metrics << std::fixed << std::setprecision(12)
            << "{\n  \"schema_version\": 2,\n  \"software\": \"lncrna_mers_v2\",\n"
            << "  \"software_version\": \"2-preview.3\",\n"
            << "  \"semantic_fasta_hash_encoding\": \"FragIncRNA-FASTA-semantic-v1\",\n"
            << "  \"compiler_version\": " << json_string(__VERSION__) << ",\n"
            << "  \"config_path\": " << json_string(o.config_path.string()) << ",\n"
            << "  \"config_sha256\": " << json_string(config_sha) << ",\n"
            << "  \"reference\": " << json_string(o.reference.string()) << ",\n"
            << "  \"reference_sequence_fnv1a64\": " << json_string(hex64(info.sequence_hash)) << ",\n"
            << "  \"reference_sequence_sha256\": " << json_string(info.sequence_sha256) << ",\n"
            << "  \"reference_repeat_bed_sha256\": " << json_string(reference_bed_sha) << ",\n"
            << "  \"queries\": " << json_string(o.queries.string()) << ",\n"
            << "  \"query_sequence_fnv1a64\": " << json_string(hex64(query_sequence_hash)) << ",\n"
            << "  \"query_sequence_sha256\": " << json_string(query_sha.hex_digest()) << ",\n"
            << "  \"query_repeat_bed_sha256\": " << json_string(query_bed_sha) << ",\n"
            << "  \"exon_blocks_sha256\": " << json_string(exon_blocks_sha) << ",\n"
            << "  \"expected_loci_bed_sha256\": " << json_string(expected_loci_sha) << ",\n"
            << "  \"threshold_table_sha256\": " << json_string(threshold_table_sha) << ",\n"
            << "  \"mode\": " << json_string(o.mode) << ",\n"
            << "  \"strand_mode\": " << json_string(o.strand) << ",\n"
            << "  \"repeat_policy\": " << json_string(o.repeat) << ",\n"
            << "  \"mask_target\": " << json_string(o.mask_target) << ",\n"
            << "  \"reference_repeat_bed\": " << json_string(o.reference_repeat_bed.string()) << ",\n"
            << "  \"query_repeat_bed\": " << json_string(o.query_repeat_bed.string()) << ",\n"
            << "  \"repeat_classes\": " << json_string(o.repeat_classes.string()) << ",\n"
            << "  \"reference_hit_classes\": " << json_string(o.reference_hit_classes.string()) << ",\n"
            << "  \"max_repeat_intervals\": " << o.max_repeat_intervals << ",\n"
            << "  \"exon_blocks\": " << json_string(o.exon_blocks.string()) << ",\n"
            << "  \"max_exon_blocks\": " << o.max_exon_blocks << ",\n"
            << "  \"discovered_loci\": " << json_string(o.loci_output.string()) << ",\n"
            << "  \"locus_gap\": " << o.locus_gap << ",\n"
            << "  \"max_locus_span\": " << o.max_locus_span << ",\n"
            << "  \"max_locus_rows\": " << o.max_locus_rows << ",\n"
            << "  \"expected_loci_bed\": " << json_string(o.expected_loci_bed.string()) << ",\n"
            << "  \"expected_loci_output\": " << json_string(o.expected_loci_output.string()) << ",\n"
            << "  \"max_expected_loci\": " << o.max_expected_loci << ",\n"
            << "  \"max_query_bases\": " << o.max_query_bases << ",\n"
            << "  \"threshold_source\": " << json_string(o.threshold_source) << ",\n"
            << "  \"threshold_table\": " << json_string(o.threshold_table.string()) << ",\n"
            << "  \"training_control_set\": " << json_string(o.training_control_set) << ",\n"
            << "  \"evaluation_control_set\": " << json_string(o.evaluation_control_set) << ",\n"
            << "  \"basic_threshold\": " << o.basic << ",\n"
            << "  \"strict_threshold\": " << o.strict << ",\n"
            << "  \"local_window\": " << o.local_window << ",\n"
            << "  \"max_position_rows_per_batch\": " << o.max_position_rows << ",\n"
            << "  \"positions\": " << json_string(o.positions.string()) << ",\n"
            << "  \"fragments\": " << json_string(o.fragments.string()) << ",\n"
            << "  \"k\": " << o.k << ",\n  \"fragment_size\": " << o.fragment_size << ",\n"
            << "  \"bin_count\": " << info.bin_count << ",\n"
            << "  \"hash_functions\": " << o.hash_functions << ",\n"
            << "  \"index_format_version\": " << (ibf ? 3 : 0) << ",\n"
            << "  \"index_provenance\": " << json_string(index_provenance) << ",\n"
            << "  \"index_path\": " << json_string(index_path.string()) << ",\n"
            << "  \"index_file_sha256\": " << json_string(index_sha) << ",\n"
            << "  \"configured_per_bin_fpr\": " << o.fpr << ",\n"
            << "  \"requested_per_bin_fpr_cap\": " << o.requested_per_bin_fpr << ",\n"
            << "  \"target_genome_wide_fpr\": ";
    if (o.mode == "ibf" || o.genome_fpr_given) metrics << o.genome_fpr; else metrics << "null";
    metrics << ",\n  \"observed_occupancy_union_estimate\": ";
    if (o.mode == "ibf") metrics << occupancy_estimate; else metrics << "null";
    metrics << ",\n  \"occupancy_guard_seconds\": ";
    if (o.mode == "ibf") metrics << occupancy_guard_seconds; else metrics << "null";
    metrics << ",\n  \"peak_rss_after_occupancy_guard_bytes\": ";
    if (o.mode == "ibf") metrics << occupancy_guard_peak; else metrics << "null";
    metrics << ",\n  \"fragment_output_seconds\": " << fragment_seconds
            << ",\n  \"peak_rss_after_fragment_output_bytes\": " << fragment_peak
            << ",\n"
            << "  \"nominal_any_bin_fpr_independent_bins\": " << nominal_any_bin_fpr << ",\n"
            << "  \"primary_match_source\": " << json_string(o.mode == "ibf" ? "ibf" : "exact") << ",\n"
            << "  \"calibration_bins\": " << json_string(o.calibration_bins.string()) << ",\n"
            << "  \"bin_bits\": " << info.bin_bits << ",\n"
            << "  \"distinct_bin_kmers_sum\": ";
    if (ibf) metrics << std::accumulate(info.bin_distinct.begin(), info.bin_distinct.end(), std::size_t{0});
    else metrics << "null";
    metrics << ",\n  \"occupied_bin_bits_sum\": ";
    if (o.mode == "calibrate") metrics << occupied_bits_total; else metrics << "null";
    metrics << ",\n  \"bit_occupancy_measurement_seconds\": ";
    if (o.mode == "calibrate") metrics << occupancy_seconds; else metrics << "null";
    metrics << ",\n  \"peak_rss_after_bit_occupancy_bytes\": ";
    if (o.mode == "calibrate") metrics << occupancy_peak; else metrics << "null";
    metrics << ",\n"
            << "  \"max_calibration_pairs_per_batch\": " << o.max_calibration_pairs << ",\n"
            << "  \"max_ibf_bin_pairs_per_batch\": " << o.max_ibf_bin_pairs << ",\n"
            << "  \"threads\": 1,\n"
            << "  \"total_reference_bases\": " << info.total_bases << ",\n"
            << "  \"valid_reference_windows\": " << info.valid_windows << ",\n"
            << "  \"query_count\": " << totals.queries << ",\n"
            << "  \"valid_query_windows\": " << totals.valid << ",\n"
            << "  \"discovered_locus_rows\": ";
    if (o.loci_output.empty()) metrics << "null"; else metrics << totals.discovered_locus_rows;
    metrics << ",\n  \"expected_locus_matched_windows\": ";
    if (o.expected_loci_bed.empty()) metrics << "null"; else metrics << totals.expected_locus_matched;
    metrics << ",\n  \"expected_locus_recovery\": ";
    if (o.expected_loci_bed.empty() || !totals.expected_valid_windows) metrics << "null";
    else metrics << static_cast<double>(totals.expected_locus_matched) / totals.expected_valid_windows;
    metrics << ",\n"
            << "  \"exon_eligible_windows\": ";
    if (o.exon_blocks.empty()) metrics << "null"; else metrics << totals.exon_eligible;
    metrics << ",\n  \"exon_matched_windows\": ";
    if (o.exon_blocks.empty()) metrics << "null"; else metrics << totals.exon_matched;
    metrics << ",\n  \"exon_aware_recovery\": ";
    if (o.exon_blocks.empty() || !totals.exon_eligible) metrics << "null";
    else metrics << static_cast<double>(totals.exon_matched) / totals.exon_eligible;
    metrics << ",\n"
            << "  \"tested_absent_windows\": " << totals.absent << ",\n"
            << "  \"exact_positive_windows\": " << totals.exact_positive << ",\n"
            << "  \"false_positive_windows\": ";
    if (o.mode == "exact" || o.mode == "ibf") metrics << "null"; else metrics << totals.fp;
    metrics << ",\n  \"false_negative_windows\": ";
    if (o.mode == "exact" || o.mode == "ibf") metrics << "null"; else metrics << totals.fn;
    metrics << ",\n  \"measured_false_negative_rate\": ";
    if (o.mode == "exact" || o.mode == "ibf" || !totals.exact_positive) metrics << "null";
    else metrics << static_cast<double>(totals.fn) / totals.exact_positive;
    metrics << ",\n  \"tested_absent_bin_pairs\": ";
    if (o.mode == "calibrate") metrics << absent_bin_pairs; else metrics << "null";
    metrics << ",\n  \"false_positive_bin_pairs\": ";
    if (o.mode == "calibrate") metrics << fp_bin_pairs; else metrics << "null";
    metrics << ",\n  \"false_negative_bin_pairs\": ";
    if (o.mode == "calibrate") metrics << fn_bin_pairs; else metrics << "null";
    metrics << ",\n  \"measured_per_bin_fpr\": ";
    if (o.mode == "calibrate" && absent_bin_pairs) metrics << static_cast<double>(fp_bin_pairs) / absent_bin_pairs;
    else metrics << "null";
    metrics << ",\n"
            << "  \"measured_genome_wide_fpr\": ";
    if (o.mode == "exact" || o.mode == "ibf" || totals.absent == 0) metrics << "null";
    else metrics << static_cast<double>(totals.fp) / totals.absent;
    metrics << ",\n  \"reference_inspection_seconds\": " << inspect_seconds
            << ",\n  \"peak_rss_after_reference_inspection_bytes\": " << inspect_peak
            << ",\n  \"index_construction_seconds\": " << index_seconds
            << ",\n  \"index_load_seconds\": " << index_load_seconds
            << ",\n  \"index_save_seconds\": " << index_save_seconds
            << ",\n  \"peak_rss_after_index_bytes\": " << index_peak
            << ",\n  \"query_seconds\": " << totals.query_seconds
            << ",\n  \"exact_scan_seconds\": " << totals.exact_seconds
            << ",\n  \"result_output_seconds\": " << totals.result_seconds
            << ",\n  \"reference_valid_window_throughput_per_second\": "
            << (index_seconds > 0 ? info.valid_windows / index_seconds : 0.0)
            << ",\n  \"query_valid_window_throughput_per_second\": "
            << (totals.query_seconds > 0 ? totals.valid / totals.query_seconds : 0.0)
            << ",\n  \"total_seconds\": " << elapsed(run_start)
            << ",\n  \"peak_resident_bytes\": " << peak_rss_bytes()
            << ",\n  \"disk_index_bytes\": " << disk_index_bytes << "\n}\n";
    std::cerr << "Wrote " << totals.queries << " query rows; elapsed=" << elapsed(run_start)
              << " s; peak RSS=" << peak_rss_bytes() << " bytes\n";
    return 0;
}
} // namespace

int main(int argc, char ** argv)
{
    try { return run(parse_options(argc, argv)); }
    catch (std::exception const & e) { std::cerr << "Error: " << e.what() << '\n'; return 1; }
}
