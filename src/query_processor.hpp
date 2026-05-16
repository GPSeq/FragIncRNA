#pragma once

#include "config.hpp"
#include "reference_index.hpp"

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
    /*
    * @fn QueryProcessor
    * @brief Creates a query processor for one reference index.
    * @signature QueryProcessor(Config const & cfg, ReferenceIndex & index, std::string const & ref_name);
    * @param cfg: application configuration.
    * @param index: reference index used to count query k-mer hits.
    * @param ref_name: reference name used in output headers and file names.
    * @throws None.
    * @return None.
    */
    QueryProcessor(Config const & cfg,
                   ReferenceIndex & index,
                   std::string const & ref_name);

    /*
    * @fn run_fill_results_col
    * @brief Processes all queries for combined output mode and fills one reference column in the result matrix.
    * @signature void run_fill_results_col(std::size_t ref_idx, std::vector<std::vector<RefResult>> & results) const;
    * @param ref_idx: reference column index to fill.
    * @param results: result matrix indexed by query and reference.
    * @throws std::runtime_error when output files cannot be opened or result dimensions do not match the query file.
    * @return None.
    */
    void run_fill_results_col(std::size_t ref_idx,
                              std::vector<std::vector<RefResult>> & results) const;

    /*
    * @fn run_write_per_ibf
    * @brief Processes all queries for per-reference output mode and streams one TSV result file.
    * @signature void run_write_per_ibf(std::filesystem::path const & out_path) const;
    * @param out_path: path to the per-reference TSV output file.
    * @throws std::runtime_error when an output file cannot be opened or sequence processing fails.
    * @return None.
    */
    void run_write_per_ibf(std::filesystem::path const & out_path) const;

private:
    Config      cfg_;
    ReferenceIndex & index_;
    std::string ref_name_;
};
