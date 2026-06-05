cmake_minimum_required(VERSION 3.20)

add_executable(negative_set
    src/main.cpp
)

target_link_libraries(negative_set
    PRIVATE
        lncrna_mers_lib
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(negative_set PRIVATE -Wall -Wextra -Wpedantic)
endif()
