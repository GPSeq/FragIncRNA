#include "bloom_sizing.hpp"
#include <cassert>
#include <iostream>

int main() {
    using namespace fraglncrna;
    // Do not compile away independent checks in Release builds.
    auto require = [](bool value) { if (!value) throw std::runtime_error("Bloom sizing assertion failed"); };
    require(union_independent(0,1) == 0);
    require(union_independent(100,0) == 0);
    require(union_independent(100,1) == 1);
    require(std::abs(union_independent(1,.01)-.01) < 1e-16);
    auto tiny = union_independent(1000000000,1e-20);
    require(std::abs(tiny-1e-11) < 1e-20);
    require(union_independent(std::numeric_limits<std::size_t>::max(),.01) == 1);
    require(union_bound(100,.01) == 1);
    require(planned_bin_bits(80000,3,.01) == 989134);
    require(planned_index_bytes(65,1024) == 16384);
    require(planned_bin_bits(0,3,.01) == 1024);
    for (int which = 0; which < 4; ++which) {
        bool rejected = false;
        try {
            if (which == 0) planned_index_bytes(std::numeric_limits<std::size_t>::max(),1024);
            if (which == 1) planned_bin_bits(80000,1,std::numeric_limits<double>::denorm_min());
            if (which == 2) union_independent(10,std::numeric_limits<double>::quiet_NaN());
            if (which == 3) planned_bin_bits(10,6,.01);
        } catch (std::runtime_error const &) { rejected = true; }
        require(rejected);
    }
    std::cout << "Bloom numeric boundary checks passed\n";
}
