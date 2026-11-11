#pragma once

/**
 * @file reference_index.hpp
 * @brief IBF/HIBF index of one reference for the legacy executable `lncrna_mers`.
 */
#include "config.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#if defined(LNCRNA_MERS_HAS_HIBF)
namespace seqan::hibf
{
class hierarchical_interleaved_bloom_filter;
}
#endif

template <typename sequence_t>
class ReferenceIndex
{
public:
    /**
     * @brief Builds an IBF or HIBF index for the fragments of one reference.
     * @param ref_name reference name associated with the index.
     * @param fragments DNA fragments to insert into the index.
     * @param cfg application configuration controlling index type and parameters.
     * @throws std::runtime_error when fragments are empty, k-mer settings are invalid, or HIBF support is unavailable.
     */
    ReferenceIndex(std::string ref_name,
                   std::vector<sequence_t> const & fragments,
                   Config const & cfg);

    /**
     * @brief Destroys the reference index and releases owned index storage.
     */
    ~ReferenceIndex();

    /**
     * @brief Returns the reference name associated with this index.
     * @return Reference name.
     */
    [[nodiscard]] std::string const & ref_name() const noexcept
    {
        return ref_name_;
    }

    /**
     * @brief Returns the number of user bins represented by the reference fragments.
     * @return Number of user bins.
     */
    [[nodiscard]] std::size_t bin_count() const noexcept;

    /**
     * @brief Counts how many index bins match each k-mer in a query sequence.
     * @param seq query DNA sequence to search against the index.
     * @throws std::runtime_error when HIBF search is requested in a build without HIBF support.
     * @return Per-k-mer hit counts.
     */
    [[nodiscard]] std::vector<std::size_t> count_query_kmer_hits(sequence_t const & seq) const;

    /**
     * @brief Returns the file suffix used when serializing the selected index type.
     * @return Index file suffix.
     */
    [[nodiscard]] std::string index_file_suffix() const;

    [[nodiscard]] std::vector<std::filesystem::path>
    index_file_paths(std::filesystem::path const & base_path) const;

    /**
     * @brief Serializes the selected index to a binary archive file.
     * @param out_path path to the binary index output file.
     * @throws std::runtime_error when the file cannot be opened or HIBF support is unavailable.
     */
    void store_to(std::filesystem::path const & out_path) const;

private:
    std::string ref_name_;
    Config cfg_;
    std::size_t user_bin_count_{0};

    std::unique_ptr<seqan3::interleaved_bloom_filter<>> ibf_;
    std::vector<std::unique_ptr<seqan3::interleaved_bloom_filter<>>> ibfs_;
#if defined(LNCRNA_MERS_HAS_HIBF)
    std::unique_ptr<seqan::hibf::hierarchical_interleaved_bloom_filter> hibf_;
#endif

    /**
     * @brief Builds a flat interleaved Bloom filter from the reference fragments.
     * @param fragments DNA fragments to insert into the IBF.
     * @throws std::runtime_error when k-mer settings are invalid.
     */
    void build_ibf(std::vector<sequence_t> const & fragments);
    void build_multi_ibf(std::vector<sequence_t> const & fragments);
#if defined(LNCRNA_MERS_HAS_HIBF)
    /**
     * @brief Builds a hierarchical interleaved Bloom filter from the reference fragments.
     * @param fragments DNA fragments to insert into the HIBF.
     * @throws std::runtime_error when HIBF configuration validation fails.
     */
    void build_hibf(std::vector<sequence_t> const & fragments);
#endif
};

using ReferenceIndexDna5 = ReferenceIndex<seqan3::dna5_vector>;
using ReferenceIndexDna4 = ReferenceIndex<seqan3::dna4_vector>;
