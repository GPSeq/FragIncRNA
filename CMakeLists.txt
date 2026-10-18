cmake_minimum_required(VERSION 3.20)

project(lncrna_mers LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
enable_testing()

set(seqan3_git_tag "3.4.1")
message(STATUS "Fetching SeqAn3 ${seqan3_git_tag}")

FetchContent_Declare(
    seqan3_fetch_content
    GIT_REPOSITORY "https://github.com/seqan/seqan3.git"
    GIT_TAG        "${seqan3_git_tag}"
)

FetchContent_MakeAvailable(seqan3_fetch_content)

option(LNCRNA_MERS_ENABLE_HIBF "Build optional HIBF support." ON)
set(HIBF_LOCAL_SOURCE_DIR "" CACHE PATH "Optional local checkout of https://github.com/lutfia95/hibf.git")

set(LNCRNA_MERS_HAS_HIBF OFF)
if (LNCRNA_MERS_ENABLE_HIBF)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12)
        message(WARNING "HIBF support disabled: GCC ${CMAKE_CXX_COMPILER_VERSION} is too old. HIBF requires GCC >= 12.")
    else()
        if (HIBF_LOCAL_SOURCE_DIR)
            message(STATUS "Using local HIBF checkout: ${HIBF_LOCAL_SOURCE_DIR}")
            FetchContent_Declare(
                hibf_fetch_content
                SOURCE_DIR "${HIBF_LOCAL_SOURCE_DIR}"
            )
        else()
            message(STATUS "Fetching HIBF from lutfia95/hibf")
            FetchContent_Declare(
                hibf_fetch_content
                GIT_REPOSITORY "https://github.com/lutfia95/hibf.git"
                GIT_TAG        "main"
            )
        endif()

        FetchContent_MakeAvailable(hibf_fetch_content)
        set(LNCRNA_MERS_HAS_HIBF ON)
    endif()
endif()

add_library(lncrna_mers_lib
    src/config.cpp
    src/logger.cpp
    src/fragmenter.cpp
    src/reference_index.cpp
    src/query_processor.cpp
)

target_include_directories(lncrna_mers_lib PUBLIC src)

target_link_libraries(lncrna_mers_lib
    PUBLIC
        seqan3::seqan3
)

if (LNCRNA_MERS_HAS_HIBF)
    target_link_libraries(lncrna_mers_lib PUBLIC seqan::hibf)
    target_compile_definitions(lncrna_mers_lib PUBLIC LNCRNA_MERS_HAS_HIBF=1)
endif()

add_executable(lncrna_mers
    src/main.cpp
)

find_package(ZLIB REQUIRED)
add_executable(lncrna_mers_v2 src/kmer_v2.cpp)
target_link_libraries(lncrna_mers_v2 PRIVATE seqan3::seqan3 ZLIB::ZLIB)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(lncrna_mers_v2 PRIVATE -Wall -Wextra -Wpedantic)
endif()

target_link_libraries(lncrna_mers
    PRIVATE
        lncrna_mers_lib
)

add_subdirectory(negative_set)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(lncrna_mers_lib PRIVATE -Wall -Wextra -Wpedantic)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(lncrna_mers_lib PRIVATE -O3)
endif()

if (BUILD_TESTING)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_executable(lncrna_mers_tests
        tests/unit/test_main.cpp
    )

    target_link_libraries(lncrna_mers_tests
        PRIVATE
            lncrna_mers_lib
    )

    add_test(NAME lncrna_mers_tests COMMAND lncrna_mers_tests)
    add_test(NAME kmer_v2_integration
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_kmer_v2.py
                     $<TARGET_FILE:lncrna_mers_v2>)
    add_test(NAME kmer_v2_toml
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_kmer_v2_toml.py
                     $<TARGET_FILE:lncrna_mers_v2>)
    add_test(NAME per_bin_ibf
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_per_bin_ibf.py
                     $<TARGET_FILE:lncrna_mers_v2>)
    add_test(NAME legacy_score_regression
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_legacy_score.py
                     $<TARGET_FILE:lncrna_mers>)
    add_test(NAME fasta_sampling
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_sample_fasta.py
                     ${CMAKE_SOURCE_DIR}/scripts/sample_fasta.py)
    add_test(NAME threshold_fitting
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_fit_kmer_thresholds.py)
    add_test(NAME multi_reference_merge
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_merge_kmer_v2.py
                     $<TARGET_FILE:lncrna_mers_v2>)
    add_test(NAME kmer_v2_provenance_loci
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_kmer_v2_provenance_loci.py
                     $<TARGET_FILE:lncrna_mers_v2>)
    add_test(NAME shared_legacy_windows
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_shared_legacy_windows.py
                     $<TARGET_FILE:lncrna_mers> $<TARGET_FILE:negative_set>)
endif()
