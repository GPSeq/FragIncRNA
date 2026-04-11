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
    ReferenceIndex(std::string ref_name,
                   std::vector<seqan3::dna5_vector> const & fragments,
                   Config const & cfg);
    ~ReferenceIndex();

    [[nodiscard]] std::string const & ref_name() const noexcept
    {
        return ref_name_;
    }

    [[nodiscard]] std::size_t bin_count() const noexcept;
    [[nodiscard]] std::vector<std::size_t> count_query_kmer_hits(seqan3::dna5_vector const & seq) const;
    [[nodiscard]] std::string index_file_suffix() const;
    void store_to(std::filesystem::path const & out_path) const;

private:
    std::string ref_name_;
    Config cfg_;
    std::size_t user_bin_count_{0};

    std::unique_ptr<seqan3::interleaved_bloom_filter<>> ibf_;
#if defined(LNCRNA_MERS_HAS_HIBF)
    std::unique_ptr<seqan::hibf::hierarchical_interleaved_bloom_filter> hibf_;
#endif

    void build_ibf(std::vector<seqan3::dna5_vector> const & fragments);
#if defined(LNCRNA_MERS_HAS_HIBF)
    void build_hibf(std::vector<seqan3::dna5_vector> const & fragments);
#endif
};
