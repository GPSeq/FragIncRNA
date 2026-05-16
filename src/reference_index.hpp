#pragma once

#include "config.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#if defined(LNCRNA_MERS_HAS_HIBF)
namespace seqan::hibf
{
class hierarchical_interleaved_bloom_filter;
}
#endif

class ReferenceIndex
{
public:
    /*
    * @fn ReferenceIndex
    * @brief Builds an IBF or HIBF index for the fragments of one reference.
    * @signature ReferenceIndex(std::string ref_name, std::vector<seqan3::dna5_vector> const & fragments, Config const & cfg);
    * @param ref_name: reference name associated with the index.
    * @param fragments: DNA fragments to insert into the index.
    * @param cfg: application configuration controlling index type and parameters.
    * @throws std::runtime_error when fragments are empty, k-mer settings are invalid, or HIBF support is unavailable.
    * @return None.
    */
    ReferenceIndex(std::string ref_name,
                   std::vector<seqan3::dna5_vector> const & fragments,
                   Config const & cfg);

    /*
    * @fn ~ReferenceIndex
    * @brief Destroys the reference index and releases owned index storage.
    * @signature ~ReferenceIndex();
    * @param None.
    * @throws None.
    * @return None.
    */
    ~ReferenceIndex();

    /*
    * @fn ref_name
    * @brief Returns the reference name associated with this index.
    * @signature std::string const & ref_name() const noexcept;
    * @param None.
    * @throws None.
    * @return Reference name.
    */
    [[nodiscard]] std::string const & ref_name() const noexcept
    {
        return ref_name_;
    }

    /*
    * @fn bin_count
    * @brief Returns the number of user bins represented by the reference fragments.
    * @signature std::size_t bin_count() const noexcept;
    * @param None.
    * @throws None.
    * @return Number of user bins.
    */
    [[nodiscard]] std::size_t bin_count() const noexcept;

    /*
    * @fn count_query_kmer_hits
    * @brief Counts how many index bins match each k-mer in a query sequence.
    * @signature std::vector<std::size_t> count_query_kmer_hits(seqan3::dna5_vector const & seq) const;
    * @param seq: query DNA sequence to search against the index.
    * @throws std::runtime_error when HIBF search is requested in a build without HIBF support.
    * @return Per-k-mer hit counts.
    */
    [[nodiscard]] std::vector<std::size_t> count_query_kmer_hits(seqan3::dna5_vector const & seq) const;

    /*
    * @fn index_file_suffix
    * @brief Returns the file suffix used when serializing the selected index type.
    * @signature std::string index_file_suffix() const;
    * @param None.
    * @throws None.
    * @return Index file suffix.
    */
    [[nodiscard]] std::string index_file_suffix() const;

    /*
    * @fn store_to
    * @brief Serializes the selected index to a binary archive file.
    * @signature void store_to(std::filesystem::path const & out_path) const;
    * @param out_path: path to the binary index output file.
    * @throws std::runtime_error when the file cannot be opened or HIBF support is unavailable.
    * @return None.
    */
    void store_to(std::filesystem::path const & out_path) const;

private:
    std::string ref_name_;
    Config cfg_;
    std::size_t user_bin_count_{0};

    std::unique_ptr<seqan3::interleaved_bloom_filter<>> ibf_;
#if defined(LNCRNA_MERS_HAS_HIBF)
    std::unique_ptr<seqan::hibf::hierarchical_interleaved_bloom_filter> hibf_;
#endif

    /*
    * @fn build_ibf
    * @brief Builds a flat interleaved Bloom filter from the reference fragments.
    * @signature void build_ibf(std::vector<seqan3::dna5_vector> const & fragments);
    * @param fragments: DNA fragments to insert into the IBF.
    * @throws std::runtime_error when k-mer settings are invalid.
    * @return None.
    */
    void build_ibf(std::vector<seqan3::dna5_vector> const & fragments);
#if defined(LNCRNA_MERS_HAS_HIBF)
    /*
    * @fn build_hibf
    * @brief Builds a hierarchical interleaved Bloom filter from the reference fragments.
    * @signature void build_hibf(std::vector<seqan3::dna5_vector> const & fragments);
    * @param fragments: DNA fragments to insert into the HIBF.
    * @throws std::runtime_error when HIBF configuration validation fails.
    * @return None.
    */
    void build_hibf(std::vector<seqan3::dna5_vector> const & fragments);
#endif
};
