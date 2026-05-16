#include "fragmenter.hpp"

#include "logger.hpp"

#include <algorithm>
#include <stdexcept>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

/*
* @fn Fragmenter
* @brief Creates a fragmenter using the supplied application configuration.
* @signature Fragmenter::Fragmenter(Config const & cfg);
* @param cfg: application configuration that controls fragment size and output behavior.
* @throws None.
* @return None.
*/
Fragmenter::Fragmenter(Config const & cfg)
    : cfg_{cfg}
{}

/*
* @fn fragment_reference
* @brief Reads one reference sequence file and returns overlapping DNA fragments.
* @signature std::vector<seqan3::dna5_vector> Fragmenter::fragment_reference(std::filesystem::path const & ref_path, std::string const & ref_id) const;
* @param ref_path: path to a FASTA or FASTA.gz reference file.
* @param ref_id: reference identifier used for logging and optional fragment output names.
* @throws std::runtime_error when fragment_size is invalid or sequence IO fails.
* @return Vector of overlapping DNA fragments.
*/
std::vector<seqan3::dna5_vector>
Fragmenter::fragment_reference(std::filesystem::path const & ref_path,
                               std::string const & ref_id) const
{
    if (cfg_.fragment_size < 4)
        throw std::runtime_error{"fragment_size must be >= 4 (to allow 3 bp overlap)."};

    std::vector<seqan3::dna5_vector> fragments;

    Logger::info("Reading reference file: " + ref_path.string());

    seqan3::sequence_file_input ref_in{ref_path};

    std::size_t frag_len = cfg_.fragment_size;
    std::size_t overlap = 3;
    std::size_t step = frag_len > overlap ? frag_len - overlap : 1;

    bool write_frags = cfg_.store_fragments;
    std::filesystem::path frag_out_path;
    if (write_frags)
    {
        std::filesystem::create_directories(cfg_.output_dir);
        frag_out_path = cfg_.output_dir / (ref_id + "_fragments.fasta");
    }

    if (write_frags)
    {
        seqan3::sequence_file_output frag_out{frag_out_path};

        std::size_t global_frag_idx = 0;
        for (auto & record : ref_in)
        {
            auto const & seq = record.sequence();
            std::size_t seq_len = seq.size();
            std::size_t pos = 0;

            while (pos < seq_len)
            {
                std::size_t remaining = seq_len - pos;
                std::size_t len = remaining >= frag_len ? frag_len : remaining;

                seqan3::dna5_vector frag;
                frag.resize(len);
                std::copy_n(seq.begin() + static_cast<std::ptrdiff_t>(pos),
                            static_cast<std::ptrdiff_t>(len),
                            frag.begin());

                fragments.push_back(frag);

                std::string frag_id = ref_id + "_frag" + std::to_string(global_frag_idx);
                frag_out.emplace_back(frag, frag_id);

                ++global_frag_idx;
                if (remaining <= step)
                    break;
                pos += step;
            }
        }
    }
    else
    {
        for (auto & record : ref_in)
        {
            auto const & seq = record.sequence();
            std::size_t seq_len = seq.size();
            std::size_t pos = 0;

            while (pos < seq_len)
            {
                std::size_t remaining = seq_len - pos;
                std::size_t len = remaining >= frag_len ? frag_len : remaining;

                seqan3::dna5_vector frag;
                frag.resize(len);
                std::copy_n(seq.begin() + static_cast<std::ptrdiff_t>(pos),
                            static_cast<std::ptrdiff_t>(len),
                            frag.begin());

                fragments.push_back(frag);

                if (remaining <= step)
                    break;
                pos += step;
            }
        }
    }

    Logger::info("Reference '" + ref_id + "' fragmented into " +
                 std::to_string(fragments.size()) + " fragments (len=" +
                 std::to_string(cfg_.fragment_size) + ", overlap=3).");

    return fragments;
}
