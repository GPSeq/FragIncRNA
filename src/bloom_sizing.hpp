#pragma once

/**
 * @file bloom_sizing.hpp
 * @brief Bloom-filter sizing and false-positive-rate formulas.
 */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace fraglncrna {
/**
 * @brief Probability that at least one of @p bins independent bins reports a false positive.
 * @details Computes 1 - (1 - p)^bins in a numerically stable way.
 * @param bins Number of Bloom-filter bins.
 * @param p Per-bin false-positive rate in [0, 1].
 * @return Any-bin false-positive probability.
 * @throws std::runtime_error if @p p is not in [0, 1].
 */
inline double union_independent(std::size_t bins, double p) {
    if (!std::isfinite(p) || p < 0 || p > 1) throw std::runtime_error("Invalid Bloom FPR");
    if (!bins || p == 0) return 0;
    if (p == 1) return 1;
    return -std::expm1(double(bins)*std::log1p(-p));
}
/**
 * @brief Union (Bonferroni) upper bound min(1, bins * p) of the any-bin false-positive rate.
 * @param bins Number of Bloom-filter bins.
 * @param p Per-bin false-positive rate in [0, 1].
 * @return Upper bound of the any-bin false-positive probability.
 * @throws std::runtime_error if @p p is not in [0, 1].
 */
inline double union_bound(std::size_t bins, double p) {
    if (!std::isfinite(p) || p < 0 || p > 1) throw std::runtime_error("Invalid Bloom FPR");
    return std::min(1.0,double(bins)*p);
}
/**
 * @brief Number of bits per bin needed for a target false-positive rate.
 * @details m = ceil(-h * n / ln(1 - fpr^(1/h))), with a minimum of 1024 bits.
 * @param load Maximum number of k-mers inserted into one bin (n; the fragment size).
 * @param hashes Number of hash functions (h), 1-5.
 * @param fpr Target per-bin false-positive rate, (0, 0.5].
 * @return Bits per bin.
 * @throws std::runtime_error on invalid settings or overflow.
 */
inline std::size_t planned_bin_bits(std::size_t load, std::size_t hashes, double fpr) {
    if (!hashes || hashes > 5 || !std::isfinite(fpr) || fpr <= 0 || fpr > .5)
        throw std::runtime_error("Invalid Bloom sizing settings");
    double bits = std::ceil(-double(hashes)*double(load)/std::log1p(-std::pow(fpr,1.0/hashes)));
    if (!std::isfinite(bits) || bits < 0 || bits > double(std::numeric_limits<std::size_t>::max()/128))
        throw std::runtime_error("IBF bit sizing overflow; increase FPR or reduce reference/bin size");
    return std::max<std::size_t>(1024,static_cast<std::size_t>(bits));
}
/**
 * @brief Size of a flat IBF in bytes; bins are padded to multiples of 64.
 * @param bins Number of real bins.
 * @param bits Bits per bin.
 * @return Index size in bytes.
 * @throws std::runtime_error on overflow.
 */
inline std::size_t planned_index_bytes(std::size_t bins, std::size_t bits) {
    if (bins > std::numeric_limits<std::size_t>::max()-63) throw std::runtime_error("IBF bin count overflow");
    auto blocks = std::max<std::size_t>(1,(bins+63)/64);
    if (bits > std::numeric_limits<std::size_t>::max()/8/blocks) throw std::runtime_error("IBF byte estimate overflow");
    return bits*blocks*8;
}
}
