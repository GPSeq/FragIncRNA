#pragma once

// SHA-256 for streaming input provenance. The algorithm follows FIPS 180-4.
#include <array>
#include <bit>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>

namespace fraglncrna
{
class Sha256
{
public:
    void update(unsigned char byte)
    {
        buffer_[used_++] = byte;
        bit_count_ += 8;
        if (used_ == buffer_.size()) { transform(); used_ = 0; }
    }

    void update(std::string_view text)
    {
        for (unsigned char byte : text) update(byte);
    }

    std::string hex_digest() const
    {
        Sha256 copy = *this;
        uint64_t original_bits = copy.bit_count_;
        copy.buffer_[copy.used_++] = 0x80;
        if (copy.used_ > 56)
        {
            while (copy.used_ < 64) copy.buffer_[copy.used_++] = 0;
            copy.transform();
            copy.used_ = 0;
        }
        while (copy.used_ < 56) copy.buffer_[copy.used_++] = 0;
        for (int i = 7; i >= 0; --i)
            copy.buffer_[copy.used_++] = static_cast<unsigned char>(original_bits >> (i * 8));
        copy.transform();
        std::ostringstream out;
        out << std::hex << std::setfill('0');
        for (uint32_t word : copy.state_) out << std::setw(8) << word;
        return out.str();
    }

private:
    static constexpr std::array<uint32_t, 64> constants_ = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
    std::array<uint32_t, 8> state_ = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                                      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    std::array<unsigned char, 64> buffer_{};
    std::size_t used_ = 0;
    uint64_t bit_count_ = 0;

    void transform()
    {
        uint32_t words[64];
        for (int i = 0; i < 16; ++i)
            words[i] = (uint32_t(buffer_[4 * i]) << 24) | (uint32_t(buffer_[4 * i + 1]) << 16) |
                       (uint32_t(buffer_[4 * i + 2]) << 8) | uint32_t(buffer_[4 * i + 3]);
        for (int i = 16; i < 64; ++i)
        {
            uint32_t x = words[i - 15], y = words[i - 2];
            words[i] = words[i - 16] + (std::rotr(x, 7) ^ std::rotr(x, 18) ^ (x >> 3)) +
                       words[i - 7] + (std::rotr(y, 17) ^ std::rotr(y, 19) ^ (y >> 10));
        }
        auto [a, b, c, d, e, f, g, h] = state_;
        for (int i = 0; i < 64; ++i)
        {
            uint32_t choose = (e & f) ^ (~e & g);
            uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
            uint32_t s1 = std::rotr(e, 6) ^ std::rotr(e, 11) ^ std::rotr(e, 25);
            uint32_t s0 = std::rotr(a, 2) ^ std::rotr(a, 13) ^ std::rotr(a, 22);
            uint32_t first = h + s1 + choose + constants_[i] + words[i];
            uint32_t second = s0 + majority;
            h = g; g = f; f = e; e = d + first;
            d = c; c = b; b = a; a = first + second;
        }
        state_[0] += a; state_[1] += b; state_[2] += c; state_[3] += d;
        state_[4] += e; state_[5] += f; state_[6] += g; state_[7] += h;
    }
};
} // namespace fraglncrna
