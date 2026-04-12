
#include "query_processor.hpp"

#include "logger.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <unordered_set>
#include <numeric>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <chrono>
#include <string>
#include <iostream>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

namespace
{

std::string sequence_kmer_to_string(auto const & seq, std::size_t pos, std::size_t kmer_size)
{
    std::string kmer;
    kmer.reserve(kmer_size);

    for (std::size_t i = 0; i < kmer_size; ++i)
        kmer.push_back(seqan3::to_char(seq[pos + i]));

    return kmer;
}

std::vector<std::string> collect_unique_matching_kmers(auto const & seq,
                                                       auto const & counts,
                                                       std::size_t kmer_size)
{
    std::vector<std::string> ibf_unique_kmers;
    std::unordered_set<std::string> seen;
    seen.reserve(counts.size());

    for (std::size_t i = 0; i < counts.size(); ++i)
    {
        if (counts[i] != 1)
            continue;

        std::string kmer = sequence_kmer_to_string(seq, i, kmer_size);
        if (seen.insert(kmer).second)
            ibf_unique_kmers.push_back(std::move(kmer));
    }

    return ibf_unique_kmers;
}

std::pair<std::string, std::string> encode_matching_kmer_hits(auto const & counts)
{
    std::ostringstream indices;
    std::ostringstream hit_counts;
    bool first = true;

    for (std::size_t i = 0; i < counts.size(); ++i)
    {
        if (counts[i] == 0)
            continue;

        if (!first)
        {
            indices << '/';
            hit_counts << '/';
        }

        indices << (i + 1);
        hit_counts << counts[i];
        first = false;
    }

    return {indices.str(), hit_counts.str()};
}

} // namespace

QueryProcessor::QueryProcessor(Config const & cfg,
                               IBFIndex & index,
                               std::string const & ref_name)
    : cfg_{cfg}
    , index_{index}
    , ref_name_{ref_name}
{}

// -------------------------------------------------------------
// Combined mode: fill one column (ref_idx) of results[q][ref_idx]
// -------------------------------------------------------------
void QueryProcessor::run_fill_results_col(std::size_t ref_idx,
                                          std::vector<std::vector<RefResult>> & results) const
{
    using clock = std::chrono::steady_clock;

    Logger::info("QueryProcessor (combined) for reference '" + ref_name_ + "'.");

    auto & ibf = index_.ibf();
    auto agent = ibf.counting_agent();

    auto unique_dir = cfg_.output_dir / "unique_mers";
    std::filesystem::create_directories(unique_dir);

    auto unique_path = unique_dir / (ref_name_ + ".tsv");
    std::ofstream unique_out(unique_path);
    if (!unique_out)
        throw std::runtime_error("Failed to open IBF unique k-mer output file: " + unique_path.string());

    auto kmer_hits_path = cfg_.output_dir / (ref_name_ + "_kmers.tsv");
    std::ofstream kmer_hits_out(kmer_hits_path);
    if (!kmer_hits_out)
        throw std::runtime_error("Failed to open k-mer hit output file: " + kmer_hits_path.string());

    unique_out << "query_index\tibf_unique_kmer\n";
    kmer_hits_out << "query_index\tkmer_indices\tkmer_counts\n";

    std::size_t total_queries = results.size();

    seqan3::sequence_file_input query_in{cfg_.query_file};

    std::size_t q = 0;
    double total_ibf_time = 0.0;

    for (auto & record : query_in)
    {
        if (q >= total_queries)
            throw std::runtime_error("More queries in file than allocated rows in results matrix.");

        auto const & seq = record.sequence();
        std::size_t total_kmers =
            seq.size() >= cfg_.kmer_size ? seq.size() - cfg_.kmer_size + 1 : 0;

        auto hash_view = seq | seqan3::views::kmer_hash(
                                    seqan3::ungapped{static_cast<uint8_t>(cfg_.kmer_size)});

        auto start  = clock::now();
        auto counts = agent.bulk_count(hash_view);
        auto stop   = clock::now();

        double dt_sec =
            std::chrono::duration<double>(stop - start).count();
        total_ibf_time += dt_sec;

        std::uint64_t match_count =
            std::accumulate(counts.begin(), counts.end(), std::uint64_t{0});
        auto ibf_unique_kmers = collect_unique_matching_kmers(seq, counts, cfg_.kmer_size);
        auto [matching_kmer_indices, matching_kmer_counts] = encode_matching_kmer_hits(counts);

        bool   pass = (match_count >= cfg_.hit_threshold);
        double pct  = (total_kmers > 0)
                      ? static_cast<double>(match_count) / total_kmers
                      : 0.0;

        if (ref_idx >= results[q].size())
            throw std::runtime_error("ref_idx out of range in results matrix.");

        results[q][ref_idx] = RefResult{match_count,
                                        static_cast<std::uint64_t>(ibf_unique_kmers.size()),
                                        pass,
                                        pct};

        if (ibf_unique_kmers.empty())
        {
            unique_out << q << '\t' << '\n';
        }
        else
        {
            for (auto const & kmer : ibf_unique_kmers)
                unique_out << q << '\t' << kmer << '\n';
        }

        kmer_hits_out << q << '\t'
                      << matching_kmer_indices << '\t'
                      << matching_kmer_counts << '\n';

        // progress
        std::ostringstream prog;
        prog << "\r[combined] ref " << (ref_idx + 1)
             << " '" << ref_name_ << "', query "
             << (q + 1) << "/" << total_queries
             << ", time=" << std::fixed << std::setprecision(3)
             << dt_sec << "s";
        Logger::print_stdout(prog.str());

        ++q;
    }

    if (q != total_queries)
        throw std::runtime_error("Fewer queries in file than rows in results matrix.");

    Logger::print_stdout("", true);
    Logger::info("Finished combined processing for '" + ref_name_ +
                 "', total IBF time: " + std::to_string(total_ibf_time) + " s.");
    Logger::info("IBF unique k-mer results written to: " + unique_path.string());
    Logger::info("K-mer hit results written to: " + kmer_hits_path.string());
}

// -------------------------------------------------------------
// Per-IBF mode: stream directly to results_<ref>.tsv
// -------------------------------------------------------------
void QueryProcessor::run_write_per_ibf(std::filesystem::path const & out_path) const
{
    using clock = std::chrono::steady_clock;

    Logger::info("QueryProcessor (per-IBF) for reference '" + ref_name_ + "'.");

    std::ofstream out(out_path);
    if (!out)
        throw std::runtime_error("Failed to open per-IBF result file: " + out_path.string());

    auto & ibf = index_.ibf();
    auto agent = ibf.counting_agent();

    auto unique_dir = cfg_.output_dir / "unique_mers";
    std::filesystem::create_directories(unique_dir);

    auto unique_path = unique_dir / (ref_name_ + ".tsv");
    std::ofstream unique_out(unique_path);
    if (!unique_out)
        throw std::runtime_error("Failed to open IBF unique k-mer output file: " + unique_path.string());

    auto kmer_hits_path = cfg_.output_dir / (ref_name_ + "_kmers.tsv");
    std::ofstream kmer_hits_out(kmer_hits_path);
    if (!kmer_hits_out)
        throw std::runtime_error("Failed to open k-mer hit output file: " + kmer_hits_path.string());

    // header
    out << "query_index"
        << '\t' << ref_name_ << "_count"
        << '\t' << ref_name_ << "_ibf_unique_kmer"
        << '\t' << ref_name_ << "_pass"
        << '\t' << ref_name_ << "_pct"
        << '\n';

    unique_out << "query_index\tibf_unique_kmer\n";
    kmer_hits_out << "query_index\tkmer_indices\tkmer_counts\n";

    seqan3::sequence_file_input query_in{cfg_.query_file};

    std::size_t q = 0;
    double total_ibf_time = 0.0;

    for (auto & record : query_in)
    {
        auto const & seq = record.sequence();

        std::size_t total_kmers =
            seq.size() >= cfg_.kmer_size ? seq.size() - cfg_.kmer_size + 1 : 0;

        auto hash_view = seq | seqan3::views::kmer_hash(
                                    seqan3::ungapped{static_cast<uint8_t>(cfg_.kmer_size)});

        auto start  = clock::now();
        auto counts = agent.bulk_count(hash_view);
        auto stop   = clock::now();

        double dt_sec =
            std::chrono::duration<double>(stop - start).count();
        total_ibf_time += dt_sec;

        std::uint64_t match_count =
            std::accumulate(counts.begin(), counts.end(), std::uint64_t{0});
        auto ibf_unique_kmers = collect_unique_matching_kmers(seq, counts, cfg_.kmer_size);
        auto [matching_kmer_indices, matching_kmer_counts] = encode_matching_kmer_hits(counts);

        bool   pass = (match_count >= cfg_.hit_threshold);
        double pct  = (total_kmers > 0)
                      ? static_cast<double>(match_count) / total_kmers
                      : 0.0;

        out << q << '\t'
            << match_count << '\t'
            << ibf_unique_kmers.size() << '\t'
            << (pass ? 1 : 0) << '\t'
            << std::fixed << std::setprecision(4) << pct
            << '\n';

        if (ibf_unique_kmers.empty())
        {
            unique_out << q << '\t' << '\n';
        }
        else
        {
            for (auto const & kmer : ibf_unique_kmers)
                unique_out << q << '\t' << kmer << '\n';
        }

        kmer_hits_out << q << '\t'
                      << matching_kmer_indices << '\t'
                      << matching_kmer_counts << '\n';

        // progress
        std::ostringstream prog;
        prog << "\r[per-IBF] ref '" << ref_name_
             << "', query " << (q + 1) 
             << ", time=" << std::fixed << std::setprecision(3)
             << dt_sec << "s";
        Logger::print_stdout(prog.str());

        ++q;
    }

    Logger::print_stdout("", true);
    Logger::info("Finished per-IBF results for '" + ref_name_ +
                 "', total IBF time: " + std::to_string(total_ibf_time) + " s.");
    Logger::info("IBF unique k-mer results written to: " + unique_path.string());
    Logger::info("K-mer hit results written to: " + kmer_hits_path.string());
}

/*
void QueryProcessor::run(std::vector<std::vector<RefResult>> & results)
{
    using clock = std::chrono::steady_clock;

    Logger::info("Starting query processing for reference '" + index_.ref_name() + "'.");

    // Open queries
    seqan3::sequence_file_input query_in{cfg_.query_file};

    auto hash_view = seqan3::views::kmer_hash(
        seqan3::ungapped{static_cast<uint8_t>(cfg_.kmer_size)});

    auto & ibf = index_.ibf();
    auto agent = ibf.counting_agent<>();

    std::size_t query_idx = 0;
    double total_ibf_time = 0.0;

    for (auto & record : query_in)
    {
        auto const & seq = record.sequence();
        std::string qid = record.id();

        std::size_t total_kmers =
            seq.size() >= cfg_.kmer_size ? seq.size() - cfg_.kmer_size + 1 : 0;

        // Optional: log into file (not stdout)
        Logger::info("Query '" + qid + "' vs '" + index_.ref_name() + "' (" +
                     std::to_string(seq.size()) + " bp, " +
                     std::to_string(total_kmers) + " k-mers).");

        auto start = clock::now();
        auto counts = agent.bulk_count(seq | hash_view);
        auto stop = clock::now();
        std::chrono::duration<double> dt = stop - start;
        double dt_sec = dt.count();
        total_ibf_time += dt_sec;

        std::uint64_t match_count =
            std::accumulate(counts.begin(), counts.end(), std::uint64_t{0});

        bool pass = match_count >= cfg_.hit_threshold;
        double pct = (total_kmers > 0)
                         ? static_cast<double>(match_count) /
                               static_cast<double>(total_kmers)
                         : 0.0;

        // store into result matrix
        if (query_idx >= results.size())
            throw std::runtime_error{"Internal error: query index out of range in results matrix."};

        if (ref_idx_ >= results[query_idx].size())
            throw std::runtime_error{"Internal error: ref index out of range in results matrix."};

        results[query_idx][ref_idx_] = RefResult{match_count, pass, pct};

        // progress bar on stdout
        std::ostringstream prog;
        prog << "\r[ref "
             << (ref_idx_ + 1) << "/" << total_refs_
             << " " << index_.ref_name()
             << "] [query "
             << (query_idx + 1) << "/" << total_queries_
             << "] last_ibf=" << std::fixed << std::setprecision(3)
             << dt_sec << "s";
        std::cout << prog.str() << std::flush;

        ++query_idx;
    }

    std::cout << std::endl;
    std::cout << "Total IBF time for " << index_.ref_name()
              << ": " << std::fixed << std::setprecision(3)
              << total_ibf_time << " s" << std::endl;

    Logger::info("Finished query processing for reference '" + index_.ref_name() + "'.");
}

*/
