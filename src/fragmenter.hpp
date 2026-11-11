#pragma once

/**
 * @file fragmenter.hpp
 * @brief Splitting of reference sequences into overlapping fragments (legacy `lncrna_mers`).
 */
#include "config.hpp"

#include <filesystem>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

template <typename sequence_t>
class Fragmenter
{
public:
    /**
     * @brief Creates a fragmenter using the supplied application configuration.
     * @param cfg application configuration that controls fragment size and output behavior.
     */
    explicit Fragmenter(Config const & cfg);

    /**
     * @brief Reads one reference sequence file and returns overlapping DNA fragments.
     * @param ref_path path to a FASTA or FASTA.gz reference file.
     * @param ref_id reference identifier used for logging and optional fragment output names.
     * @throws std::runtime_error when fragment_size is invalid or sequence IO fails.
     * @return Vector of overlapping DNA fragments.
     */
    std::vector<sequence_t>
    fragment_reference(std::filesystem::path const & ref_path,
                       std::string const & ref_id) const;

private:
    Config cfg_;
};

using FragmenterDna5 = Fragmenter<seqan3::dna5_vector>;
using FragmenterDna4 = Fragmenter<seqan3::dna4_vector>;
