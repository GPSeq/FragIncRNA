#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>

struct Config
{
    std::filesystem::path ref_dir;
    std::filesystem::path query_file;

    std::filesystem::path output_dir{"."};
    std::filesystem::path output_file{"lncrna_mers_results.tsv"};
    std::filesystem::path log_file{"lncrna_mers.log"};

    std::size_t fragment_size{};       // fragment length
    std::size_t kmer_size{14};        // k-mer length (default 14, user-defined)
    std::size_t hash_functions{3};    // default 3, user-defined

    double fpr{0.01};                 // target FPR per bin, default 0.01, user-defined

    std::uint64_t hit_threshold{};    // absolute k-mer hit threshold

    bool store_fragments{false};      // write fragments to FASTA
    bool store_ibf{false};            // serialize IBFs to disk

    bool cleanup_ibf{false};

    // true = one big TSV; false = per-IBF results_<ref>.tsv (no results in RAM)
    bool single_results_writer{true};

    std::size_t threads{0};          // 0 = auto-detect, otherwise worker count
};

Config load_config_from_toml(std::filesystem::path const & config_path);
