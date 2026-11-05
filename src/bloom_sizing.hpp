#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace fraglncrna {
inline double union_independent(std::size_t bins, double p) {
    if (!std::isfinite(p) || p < 0 || p > 1) throw std::runtime_error("Invalid Bloom FPR");
    if (!bins || p == 0) return 0;
    if (p == 1) return 1;
    return -std::expm1(double(bins)*std::log1p(-p));
}
inline double union_bound(std::size_t bins, double p) {
    if (!std::isfinite(p) || p < 0 || p > 1) throw std::runtime_error("Invalid Bloom FPR");
    return std::min(1.0,double(bins)*p);
}
inline std::size_t planned_bin_bits(std::size_t load, std::size_t hashes, double fpr) {
    if (!hashes || hashes > 5 || !std::isfinite(fpr) || fpr <= 0 || fpr > .5)
        throw std::runtime_error("Invalid Bloom sizing settings");
    double bits = std::ceil(-double(hashes)*double(load)/std::log1p(-std::pow(fpr,1.0/hashes)));
    if (!std::isfinite(bits) || bits < 0 || bits > double(std::numeric_limits<std::size_t>::max()/128))
        throw std::runtime_error("IBF bit sizing overflow; increase FPR or reduce reference/bin size");
    return std::max<std::size_t>(1024,static_cast<std::size_t>(bits));
}
inline std::size_t planned_index_bytes(std::size_t bins, std::size_t bits) {
    if (bins > std::numeric_limits<std::size_t>::max()-63) throw std::runtime_error("IBF bin count overflow");
    auto blocks = std::max<std::size_t>(1,(bins+63)/64);
    if (bits > std::numeric_limits<std::size_t>::max()/8/blocks) throw std::runtime_error("IBF byte estimate overflow");
    return bits*blocks*8;
}
}
