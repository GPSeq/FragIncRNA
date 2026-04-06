cmake_minimum_required(VERSION 3.20)

project(lncrna_mers LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
enable_testing()

set(seqan3_git_tag "3.0.3")
message(STATUS "Fetching SeqAn3 ${seqan3_git_tag}")

FetchContent_Declare(
    seqan3_fetch_content
    GIT_REPOSITORY "https://github.com/seqan/seqan3.git"
    GIT_TAG        "${seqan3_git_tag}"
)

FetchContent_MakeAvailable(seqan3_fetch_content)

add_library(lncrna_mers_lib
    src/config.cpp
    src/logger.cpp
    src/fragmenter.cpp
    src/ibf_index.cpp
    src/query_processor.cpp
)

target_include_directories(lncrna_mers_lib PUBLIC src)

target_link_libraries(lncrna_mers_lib
    PUBLIC
        seqan3::seqan3
)

add_executable(lncrna_mers
    src/main.cpp
)

target_link_libraries(lncrna_mers
    PRIVATE
        lncrna_mers_lib
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(lncrna_mers_lib PRIVATE -Wall -Wextra -Wpedantic)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(lncrna_mers_lib PRIVATE -O3)
endif()

if (BUILD_TESTING)
    add_executable(lncrna_mers_tests
        tests/unit/test_main.cpp
    )

    target_link_libraries(lncrna_mers_tests
        PRIVATE
            lncrna_mers_lib
    )

    add_test(NAME lncrna_mers_tests COMMAND lncrna_mers_tests)
endif()
