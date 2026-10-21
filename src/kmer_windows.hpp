#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <optional>
#include <stdexcept>

namespace fraglncrna
{

inline int dna_code(char c) noexcept
{
    switch (c)
    {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

struct DnaWindow
{
    std::uint64_t forward = 0;
    std::uint64_t reverse = 0;
    bool valid = false;
    [[nodiscard]] std::uint64_t canonical() const noexcept { return std::min(forward, reverse); }
};

// Coordinates are preserved: an ambiguous base invalidates every overlapping window.
class RollingDnaWindow
{
public:
    explicit RollingDnaWindow(std::size_t k) : k_{k}, mask_{mask_for(k)} {}

    void reset() { invalid_.clear(); bad_ = 0; forward_ = reverse_ = 0; }

    std::optional<DnaWindow> push(char c)
    {
        int code = dna_code(c);
        if (invalid_.size() == k_)
        {
            bad_ -= invalid_.front();
            invalid_.pop_front();
        }
        invalid_.push_back(code < 0);
        bad_ += code < 0;
        auto value = static_cast<unsigned>(std::max(0, code));
        forward_ = ((forward_ << 2) | value) & mask_;
        reverse_ = (reverse_ >> 2) | (std::uint64_t{3u - value} << (2 * (k_ - 1)));
        if (invalid_.size() < k_) return std::nullopt;
        return DnaWindow{forward_, reverse_, bad_ == 0};
    }

private:
    static std::uint64_t mask_for(std::size_t k)
    {
        if (k < 1 || k > 32) throw std::invalid_argument("k must be in [1,32]");
        return k == 32 ? ~std::uint64_t{0} : (std::uint64_t{1} << (2 * k)) - 1;
    }
    std::size_t k_;
    std::uint64_t mask_, forward_ = 0, reverse_ = 0;
    std::deque<bool> invalid_;
    std::size_t bad_ = 0;
};

} // namespace fraglncrna
