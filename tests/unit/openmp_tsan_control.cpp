// Diagnostic control, deliberately not a CTest: some GCC/libgomp combinations
// do not expose OpenMP fork/join synchronization to ThreadSanitizer.
#include <array>
#include <iostream>
#include <thread>
#include <vector>
#include <omp.h>

int main(int argc, char **) {
    std::array<int,4> values{};
    int actual = 0;
    if (argc > 1) {
        std::vector<std::thread> threads;
        for (int i = 0; i < 4; ++i) threads.emplace_back([&,i] { values[i] = i+1; });
        for (auto & thread : threads) thread.join();
        actual = 4;
    } else {
#pragma omp parallel num_threads(4)
        {
#pragma omp single
            actual = omp_get_num_threads();
#pragma omp for
            for (int i = 0; i < 4; ++i) values[i] = i+1;
        } // required implicit OpenMP barrier, no concurrent main access
    }
    for (int i = 0; i < 4; ++i) if (values[i] != i+1) return 1;
    std::cout << "synchronized control completed; workers=" << actual << '\n';
}
