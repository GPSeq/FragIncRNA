#pragma once

/**
 * @file query_processor.hpp
 * @brief Query k-mer counting and result writing of the legacy executable `lncrna_mers`.
 */
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

template <typename sequence_t>
class QueryProcessor
{
public:
    /**
     * @brief Creates a query processor for one reference index.
     * @param cfg application configuration.
     * @param index reference index used to count query k-mer hits.
     * @param ref_name reference name used in output headers and file names.
     */
    QueryProcessor(Config const & cfg,
                   ReferenceIndex<sequence_t> & index,
                   std::string const & ref_name);

    /**
     * @brief Processes all queries for combined output mode and fills one reference column in the result matrix.
     * @param ref_idx reference column index to fill.
     * @param results result matrix indexed by query and reference.
     * @throws std::runtime_error when output files cannot be opened or result dimensions do not match the query file.
     */
    void run_fill_results_col(std::size_t ref_idx,
                              std::vector<std::vector<RefResult>> & results) const;

    /**
     * @brief Processes all queries for per-reference output mode and streams one TSV result file.
     * @param out_path path to the per-reference TSV output file.
     * @throws std::runtime_error when an output file cannot be opened or sequence processing fails.
     */
    void run_write_per_ibf(std::filesystem::path const & out_path) const;

private:
    Config      cfg_;
    ReferenceIndex<sequence_t> & index_;
    std::string ref_name_;
};

using QueryProcessorDna5 = QueryProcessor<seqan3::dna5_vector>;
using QueryProcessorDna4 = QueryProcessor<seqan3::dna4_vector>;
