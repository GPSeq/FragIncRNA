#pragma once

#include "config.hpp"
#include "ibf_index.hpp"

#include <cstdint>
#include <vector>

struct RefResult
{
    std::uint64_t count{};
    std::uint64_t ibf_unique_kmers{};
    bool          pass{};
    double        pct{};
};

class QueryProcessor
{
public:
    QueryProcessor(Config const & cfg,
                   IBFIndex & index,
                   std::string const & ref_name);

    // For combined mode: fill one reference column in the results matrix.
    // results.size() == #queries, results[q].size() == #refs
    void run_fill_results_col(std::size_t ref_idx,
                              std::vector<std::vector<RefResult>> & results) const;

    // For per-IBF mode: write one TSV file results_<ref>.tsv, streaming, no big matrix.
    void run_write_per_ibf(std::filesystem::path const & out_path) const;

private:
    Config      cfg_;
    IBFIndex  & index_;
    std::string ref_name_;
};
