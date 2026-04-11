#include "reference_index.hpp"

#include "logger.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

#include <cereal/archives/binary.hpp>

#include <seqan3/search/views/kmer_hash.hpp>
#if defined(LNCRNA_MERS_HAS_HIBF)
#include <hibf/config.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>
#endif

namespace
{

std::size_t compute_flat_bin_bits(std::vector<seqan3::dna5_vector> const & fragments, Config const & cfg, std::string const & ref_name)
{
    std::size_t max_len = 0;
    for (auto const & f : fragments)
        max_len = std::max<std::size_t>(max_len, f.size());

    if (cfg.kmer_size == 0 || cfg.kmer_size > max_len)
        throw std::runtime_error{
            "kmer_size must be > 0 and <= maximum fragment length for reference '" + ref_name + "'."};

    double n = static_cast<double>(max_len >= cfg.kmer_size ? max_len - cfg.kmer_size + 1 : 1);
    double k = static_cast<double>(cfg.ibf.hash_functions);
    double p = cfg.ibf.fpr;

    double p_root = std::pow(p, 1.0 / k);
    double denom = std::log(1.0 - p_root);
    if (!std::isfinite(denom) || denom == 0.0)
        denom = -1.0;

    auto const bits = static_cast<std::size_t>(std::ceil(-k * n / denom));
    return std::max<std::size_t>(bits, 1024u);
}

std::vector<std::size_t> count_hits_with_flat_ibf(seqan3::interleaved_bloom_filter<> const & ibf,
                                                  seqan3::dna5_vector const & seq,
                                                  std::size_t const kmer_size)
{
    std::vector<std::size_t> counts;

    if (seq.size() < kmer_size)
        return counts;

    counts.reserve(seq.size() - kmer_size + 1);

    auto agent = ibf.membership_agent();
    auto hash_view = seq | seqan3::views::kmer_hash(seqan3::ungapped{static_cast<uint8_t>(kmer_size)});

    for (auto const hash : hash_view)
    {
        auto const & matches = agent.bulk_contains(hash);
        counts.push_back(std::ranges::count(matches, 1));
    }

    return counts;
}

#if defined(LNCRNA_MERS_HAS_HIBF)
std::vector<std::size_t> count_hits_with_hibf(seqan::hibf::hierarchical_interleaved_bloom_filter const & hibf,
                                              seqan3::dna5_vector const & seq,
                                              std::size_t const kmer_size)
{
    std::vector<std::size_t> counts;

    if (seq.size() < kmer_size)
        return counts;

    counts.reserve(seq.size() - kmer_size + 1);

    auto agent = hibf.membership_agent();
    auto hash_view = seq | seqan3::views::kmer_hash(seqan3::ungapped{static_cast<uint8_t>(kmer_size)});

    for (auto const hash : hash_view)
    {
        std::array<uint64_t, 1> const single_hash{static_cast<uint64_t>(hash)};
        auto const & matches = agent.membership_for(single_hash, 1u);
        counts.push_back(matches.size());
    }

    return counts;
}
#endif

} // namespace

ReferenceIndex::ReferenceIndex(std::string ref_name,
                               std::vector<seqan3::dna5_vector> const & fragments,
                               Config const & cfg)
    : ref_name_{std::move(ref_name)}
    , cfg_{cfg}
    , user_bin_count_{fragments.size()}
{
    if (fragments.empty())
        throw std::runtime_error{"No fragments supplied for reference '" + ref_name_ + "'."};

    if (cfg_.index_method == IndexMethod::ibf)
    {
        build_ibf(fragments);
    }
#if defined(LNCRNA_MERS_HAS_HIBF)
    else
    {
        build_hibf(fragments);
    }
#else
    else
    {
        throw std::runtime_error("This build does not include HIBF support. Reconfigure with a GCC >= 12 toolchain.");
    }
#endif
}

ReferenceIndex::~ReferenceIndex() = default;

std::size_t ReferenceIndex::bin_count() const noexcept
{
    return user_bin_count_;
}

std::vector<std::size_t> ReferenceIndex::count_query_kmer_hits(seqan3::dna5_vector const & seq) const
{
    if (cfg_.index_method == IndexMethod::ibf)
        return count_hits_with_flat_ibf(*ibf_, seq, cfg_.kmer_size);

#if defined(LNCRNA_MERS_HAS_HIBF)
    return count_hits_with_hibf(*hibf_, seq, cfg_.kmer_size);
#else
    throw std::runtime_error("This build does not include HIBF support.");
#endif
}

std::string ReferenceIndex::index_file_suffix() const
{
    return cfg_.index_method == IndexMethod::ibf ? ".ibf" : ".hibf";
}

void ReferenceIndex::store_to(std::filesystem::path const & out_path) const
{
    std::ofstream os(out_path, std::ios::binary);
    if (!os)
        throw std::runtime_error{"Failed to open index output file: " + out_path.string()};

    cereal::BinaryOutputArchive archive(os);

    if (cfg_.index_method == IndexMethod::ibf)
        archive(*ibf_);
#if defined(LNCRNA_MERS_HAS_HIBF)
    else
        archive(*hibf_);
#else
    else
        throw std::runtime_error("This build does not include HIBF support.");
#endif
}

void ReferenceIndex::build_ibf(std::vector<seqan3::dna5_vector> const & fragments)
{
    using seqan3::bin_count;
    using seqan3::bin_index;
    using seqan3::bin_size;
    using seqan3::hash_function_count;

    std::size_t max_len = 0;
    for (auto const & f : fragments)
        max_len = std::max<std::size_t>(max_len, f.size());

    if (cfg_.kmer_size == 0 || cfg_.kmer_size > max_len)
        throw std::runtime_error{
            "kmer_size must be > 0 and <= maximum fragment length for reference '" + ref_name_ + "'."};

    auto const bin_bits = compute_flat_bin_bits(fragments, cfg_, ref_name_);

    ibf_ = std::make_unique<seqan3::interleaved_bloom_filter<>>(
        bin_count{fragments.size()},
        bin_size{bin_bits},
        hash_function_count{cfg_.ibf.hash_functions}
    );

    auto hash_view = seqan3::views::kmer_hash(seqan3::ungapped{static_cast<uint8_t>(cfg_.kmer_size)});

    std::size_t bin_idx = 0;
    for (auto const & fragment : fragments)
    {
        for (auto const hash : fragment | hash_view)
            ibf_->emplace(static_cast<uint64_t>(hash), bin_index{bin_idx});
        ++bin_idx;
    }

    Logger::print_stdout("Built IBF for reference '" + ref_name_ + "'", true);
    Logger::info("Built IBF for reference '" + ref_name_ + "' (" +
                 std::to_string(fragments.size()) + " bins, " +
                 std::to_string(bin_bits) + " bits per bin, " +
                 std::to_string(cfg_.ibf.hash_functions) + " hash functions, " +
                 "target FPR=" + std::to_string(cfg_.ibf.fpr) + ").");
}

#if defined(LNCRNA_MERS_HAS_HIBF)
void ReferenceIndex::build_hibf(std::vector<seqan3::dna5_vector> const & fragments)
{
    auto hash_view = seqan3::views::kmer_hash(seqan3::ungapped{static_cast<uint8_t>(cfg_.kmer_size)});

    auto input_fn = [&](std::size_t const user_bin_id, seqan::hibf::insert_iterator && it)
    {
        for (auto const hash : fragments[user_bin_id] | hash_view)
            it = static_cast<uint64_t>(hash);
    };

    seqan::hibf::config hibf_cfg{};
    hibf_cfg.input_fn = input_fn;
    hibf_cfg.number_of_user_bins = fragments.size();
    hibf_cfg.number_of_hash_functions = cfg_.hibf.hash_functions;
    hibf_cfg.maximum_fpr = cfg_.hibf.maximum_fpr;
    hibf_cfg.relaxed_fpr = cfg_.hibf.relaxed_fpr;
    hibf_cfg.threads = cfg_.hibf.threads == 0 ? cfg_.threads : cfg_.hibf.threads;
    hibf_cfg.sketch_bits = static_cast<uint8_t>(cfg_.hibf.sketch_bits);
    hibf_cfg.tmax = cfg_.hibf.tmax;
    hibf_cfg.empty_bin_fraction = cfg_.hibf.empty_bin_fraction;
    hibf_cfg.alpha = cfg_.hibf.alpha;
    hibf_cfg.max_rearrangement_ratio = cfg_.hibf.max_rearrangement_ratio;
    hibf_cfg.disable_estimate_union = cfg_.hibf.disable_estimate_union;
    hibf_cfg.disable_rearrangement = cfg_.hibf.disable_rearrangement;
    hibf_cfg.validate_and_set_defaults();

    hibf_ = std::make_unique<seqan::hibf::hierarchical_interleaved_bloom_filter>(hibf_cfg);

    Logger::print_stdout("Built HIBF for reference '" + ref_name_ + "'", true);
    Logger::info("Built HIBF for reference '" + ref_name_ + "' (" +
                 std::to_string(fragments.size()) + " user bins, " +
                 std::to_string(hibf_cfg.number_of_hash_functions) + " hash functions, " +
                 "maximum_fpr=" + std::to_string(hibf_cfg.maximum_fpr) +
                 ", relaxed_fpr=" + std::to_string(hibf_cfg.relaxed_fpr) +
                 ", threads=" + std::to_string(hibf_cfg.threads) + ").");
}
#endif
