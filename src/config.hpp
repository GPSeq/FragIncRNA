#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>

enum class IndexMethod
{
    ibf,
    hibf
};

struct IBFConfig
{
    std::size_t hash_functions{3};
    double fpr{0.01};
};

struct HIBFConfig
{
    std::size_t hash_functions{2};
    double maximum_fpr{0.05};
    double relaxed_fpr{0.3};
    std::size_t threads{0}; // 0 = use Config::threads
    std::size_t sketch_bits{12};
    std::size_t tmax{0};
    double empty_bin_fraction{0.0};
    double alpha{1.2};
    double max_rearrangement_ratio{0.5};
    bool disable_estimate_union{false};
    bool disable_rearrangement{false};
};

struct Config
{
    std::filesystem::path ref_dir;
    std::filesystem::path query_file;

    std::filesystem::path output_dir{"."};
    std::filesystem::path output_file{"lncrna_mers_results.tsv"};
    std::filesystem::path log_file{"lncrna_mers.log"};

    IndexMethod index_method{IndexMethod::ibf};

    std::size_t fragment_size{};       // fragment length
    std::size_t kmer_size{14};        // k-mer length (default 14, user-defined)

    std::uint64_t hit_threshold{};    // absolute k-mer hit threshold

    bool store_fragments{false};      // write fragments to FASTA
    bool store_index{false};          // serialize the selected index to disk
    bool cleanup_index{false};        // delete serialized index after processing

    // true = one big TSV; false = per-IBF results_<ref>.tsv (no results in RAM)
    bool single_results_writer{true};

    std::size_t threads{0};          // 0 = auto-detect, otherwise worker count

    IBFConfig ibf{};
    HIBFConfig hibf{};
};

Config load_config_from_toml(std::filesystem::path const & config_path);
