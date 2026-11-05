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
option(LNCRNA_MERS_ENABLE_OPENMP "Enable v2 query/interval OpenMP parallelism" ON)
if (LNCRNA_MERS_ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS CXX)
    if (OpenMP_CXX_FOUND)
        target_link_libraries(lncrna_mers_v2 PRIVATE OpenMP::OpenMP_CXX)
    else()
        message(WARNING "OpenMP unavailable: v2 will report and use serial fallback")
    endif()
endif()
execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                OUTPUT_VARIABLE FRAGLNCRNA_GIT_SHA OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
target_compile_definitions(lncrna_mers_v2 PRIVATE FRAGLNCRNA_GIT_SHA="${FRAGLNCRNA_GIT_SHA}")
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
    add_executable(bloom_sizing_tests tests/unit/test_bloom_sizing.cpp)
    target_include_directories(bloom_sizing_tests PRIVATE src)
    add_test(NAME bloom_sizing_tests COMMAND bloom_sizing_tests)
    add_test(NAME kmer_v2_integration
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_kmer_v2.py
                     $<TARGET_FILE:lncrna_mers_v2>)
    add_test(NAME kmer_v2_toml
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_kmer_v2_toml.py
                     $<TARGET_FILE:lncrna_mers_v2>)
    option(LNCRNA_MERS_REAL_DATA_TESTS "Run optional bounded real-data shard comparisons" OFF)
    set(LNCRNA_MERS_REAL_REFERENCE "${CMAKE_SOURCE_DIR}/rerun_test/genomes/T2TMMU8v20_GCF_049350105.2_T2T-MMU8v2.0_genomic.fna.gz" CACHE FILEPATH "One real FASTA reference")
    set(LNCRNA_MERS_REAL_QUERIES "${CMAKE_SOURCE_DIR}/rerun_test/gencode.v49.lncRNA_transcripts.fa" CACHE FILEPATH "GENCODE query FASTA")
    if (LNCRNA_MERS_REAL_DATA_TESTS)
        add_test(NAME sharded_ibf_real_data
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_sharded_real_data.py
                    --program $<TARGET_FILE:lncrna_mers_v2>
                    --reference ${LNCRNA_MERS_REAL_REFERENCE} --queries ${LNCRNA_MERS_REAL_QUERIES}
                    --work-dir ${CMAKE_BINARY_DIR}/real-data-shards)
        set_tests_properties(sharded_ibf_real_data PROPERTIES LABELS "real_data;optional" SKIP_RETURN_CODE 77)
    endif()
    add_executable(trusted_archive_tests tests/unit/test_trusted_archive.cpp)
    target_include_directories(trusted_archive_tests PRIVATE src)
    add_test(NAME trusted_archive_integrity COMMAND trusted_archive_tests)
    add_test(NAME production_extensions
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_production_extensions.py
                     $<TARGET_FILE:lncrna_mers_v2>)
    add_test(NAME sharded_ibf
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_sharded_ibf.py
                     $<TARGET_FILE:lncrna_mers_v2>)
    add_test(NAME per_bin_ibf
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_per_bin_ibf.py
                     $<TARGET_FILE:lncrna_mers_v2>)
    add_test(NAME per_bin_production
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_per_bin_production.py
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
    add_test(NAME production_reference_merge
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_merge_production.py
                     $<TARGET_FILE:lncrna_mers_v2>)
    add_test(NAME production_benchmark_runner
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_production_runner.py)
    add_test(NAME command_resource_measurement
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_measure_command.py)
    add_test(NAME output_documentation_contract
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_output_documentation.py
                     $<TARGET_FILE:lncrna_mers_v2>)
    add_test(NAME kmer_v2_provenance_loci
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_kmer_v2_provenance_loci.py
                     $<TARGET_FILE:lncrna_mers_v2>)
    add_test(NAME shared_legacy_windows
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_shared_legacy_windows.py
                     $<TARGET_FILE:lncrna_mers> $<TARGET_FILE:negative_set>)
endif()
