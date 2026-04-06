#include "config.hpp"
#include "fragmenter.hpp"
#include "ibf_index.hpp"
#include "logger.hpp"
#include "query_processor.hpp"

#include <deque>
#include <exception>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include <seqan3/io/sequence_file/input.hpp>

namespace fs = std::filesystem;

static bool is_fasta_like(fs::path const & p)
{
    std::string s = p.string();
    auto ends_with = [&s](std::string_view suf)
    {
        if (s.size() < suf.size())
            return false;
        return std::string_view{s}.substr(s.size() - suf.size()) == suf;
    };

    return ends_with(".fa") || ends_with(".fna") || ends_with(".fasta") ||
           ends_with(".fa.gz") || ends_with(".fna.gz") || ends_with(".fasta.gz");
}

int main(int argc, char ** argv)
{
    fs::path config_path = "config.toml";

    try
    {
        if (argc > 2)
            throw std::runtime_error("Usage: ./lncrna_mers [config.toml]");

        if (argc == 2)
            config_path = argv[1];
    }
    catch (std::exception const & ext)
    {
        std::cerr << ext.what() << '\n';
        return 1;
    }

    Config cfg;
    try
    {
        cfg = load_config_from_toml(config_path);
        fs::create_directories(cfg.output_dir);
        auto log_path = cfg.output_dir / cfg.log_file;
        Logger::init(log_path.string());
        Logger::info("Starting lncrna_mers.");
        Logger::info("Loaded config from: " + config_path.string());

        if (!fs::exists(cfg.ref_dir) || !fs::is_directory(cfg.ref_dir))
            throw std::runtime_error("Reference directory does not exist or is not a directory: " +
                                     cfg.ref_dir.string());

        if (cfg.kmer_size > cfg.fragment_size)
            throw std::runtime_error("kmer-size must be <= fragment-size.");

        if (cfg.fpr <= 0.0 || cfg.fpr >= 1.0)
            throw std::runtime_error("FPR must be in (0,1).");

        if (cfg.threads == 0)
            cfg.threads = std::max<std::size_t>(1u, std::thread::hardware_concurrency());

        // Collect reference files
        std::vector<fs::path> ref_files;
        for (auto const & entry : fs::directory_iterator{cfg.ref_dir})
        {
            if (!entry.is_regular_file())
                continue;
            if (!is_fasta_like(entry.path()))
                continue;
            ref_files.push_back(entry.path());
        }

        if (ref_files.empty())
            throw std::runtime_error("No reference FASTA / FASTA.gz files found in: " +
                                     cfg.ref_dir.string());

        Logger::info("Found " + std::to_string(ref_files.size()) +
                     " reference files to index.");
        Logger::info("Using " + std::to_string(cfg.threads) + " worker thread(s).");

        // Pre-scan queries to count them (for combined mode)
        std::size_t total_queries = 0;
        {
            seqan3::sequence_file_input query_in{cfg.query_file};
            for (auto & rec : query_in)
                ++total_queries;
        }
        if (total_queries == 0)
            throw std::runtime_error("No queries found in file: " + cfg.query_file.string());

        Logger::info("Total queries: " + std::to_string(total_queries));

        // build ref IDs
        std::vector<std::string> ref_ids;
        ref_ids.reserve(ref_files.size());
        for (auto const & ref_path : ref_files)
        {
            auto stem      = ref_path.stem().string();
            auto inner     = fs::path{stem}.stem().string();
            std::string id = inner.empty() ? stem : inner;
            ref_ids.push_back(id);
        }

        // For combined mode, allocate results matrix
        std::vector<std::vector<RefResult>> results;
        if (cfg.single_results_writer)
        {
            results.assign(total_queries,
                           std::vector<RefResult>(ref_files.size()));
        }

        auto process_reference = [&](std::size_t r)
        {
            auto const & ref_path = ref_files[r];
            auto const & ref_id   = ref_ids[r];

            Logger::info("Processing reference: " + ref_path.string() +
                         " (id='" + ref_id + "')");

            Fragmenter fragmenter{cfg};
            auto fragments = fragmenter.fragment_reference(ref_path, ref_id);

            IBFIndex ibf_idx{ref_id, fragments, cfg};

            QueryProcessor qp{cfg, ibf_idx, ref_id};

            if (cfg.single_results_writer)
            {
                // fill column r of results
                qp.run_fill_results_col(r, results);
            }
            else
            {
                // write per-IBF TSV immediately
                auto out_path = cfg.output_dir / ("results_" + ref_id + ".tsv");
                qp.run_write_per_ibf(out_path);
            }

            // Delete IBF on disk immediately after use, if requested
            if (cfg.cleanup_ibf && cfg.store_ibf)
            {
                auto ibf_path = cfg.output_dir / (ref_id + ".ibf");
                std::error_code ec;
                if (fs::exists(ibf_path, ec))
                {
                    fs::remove(ibf_path, ec);
                    if (!ec)
                        Logger::info("Removed IBF file: " + ibf_path.string());
                    else
                        Logger::warn("Failed to remove IBF file: " + ibf_path.string() +
                                     " (" + ec.message() + ")");
                }
            }
        };

        std::size_t worker_count = std::min<std::size_t>(cfg.threads, ref_files.size());
        std::deque<std::future<void>> active_workers;

        for (std::size_t r = 0; r < ref_files.size(); ++r)
        {
            active_workers.push_back(std::async(std::launch::async, process_reference, r));
            if (active_workers.size() >= worker_count)
            {
                active_workers.front().get();
                active_workers.pop_front();
            }
        }

        while (!active_workers.empty())
        {
            active_workers.front().get();
            active_workers.pop_front();
        }

        // Combined mode: write final TSV once from results matrix
        if (cfg.single_results_writer)
        {
            auto out_path = cfg.output_dir / cfg.output_file;
            std::ofstream out(out_path);
            if (!out)
                throw std::runtime_error("Failed to open result output file: " + out_path.string());

            // header
            out << "query_index";
            for (std::size_t r = 0; r < ref_ids.size(); ++r)
            {
                out << '\t' << ref_ids[r] << "_count"
                    << '\t' << ref_ids[r] << "_ibf_unique_kmer"
                    << '\t' << ref_ids[r] << "_pass"
                    << '\t' << ref_ids[r] << "_pct";
            }
            out << '\n';

            // rows
            for (std::size_t q = 0; q < total_queries; ++q)
            {
                out << q;
                for (std::size_t r = 0; r < ref_ids.size(); ++r)
                {
                    auto const & rr = results[q][r];
                    out << '\t' << rr.count
                        << '\t' << rr.ibf_unique_kmers
                        << '\t' << (rr.pass ? 1 : 0)
                        << '\t' << std::fixed << std::setprecision(4) << rr.pct;
                }
                out << '\n';
            }

            Logger::info("Combined results written to: " + out_path.string());
        }

        Logger::info("lncrna_mers finished successfully.");
    }
    catch (std::exception const & e)
    {
        Logger::error(std::string{"Fatal error: "} + e.what());
        std::cerr << "Fatal error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
