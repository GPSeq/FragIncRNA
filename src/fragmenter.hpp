#pragma once

#include "config.hpp"

#include <filesystem>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>

class Fragmenter
{
public:
    /*
    * @fn Fragmenter
    * @brief Creates a fragmenter using the supplied application configuration.
    * @signature explicit Fragmenter(Config const & cfg);
    * @param cfg: application configuration that controls fragment size and output behavior.
    * @throws None.
    * @return None.
    */
    explicit Fragmenter(Config const & cfg);

    /*
    * @fn fragment_reference
    * @brief Reads one reference sequence file and returns overlapping DNA fragments.
    * @signature std::vector<seqan3::dna5_vector> fragment_reference(std::filesystem::path const & ref_path, std::string const & ref_id) const;
    * @param ref_path: path to a FASTA or FASTA.gz reference file.
    * @param ref_id: reference identifier used for logging and optional fragment output names.
    * @throws std::runtime_error when fragment_size is invalid or sequence IO fails.
    * @return Vector of overlapping DNA fragments.
    */
    std::vector<seqan3::dna5_vector>
    fragment_reference(std::filesystem::path const & ref_path,
                       std::string const & ref_id) const;

private:
    Config cfg_;
};
