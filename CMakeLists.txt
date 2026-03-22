cmake_minimum_required(VERSION 3.14)
project(clnj_pipeline CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

# ── Eigen (header-only linear algebra) ────────────────────────────
FetchContent_Declare(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    URL_HASH MD5=4c527a9171d71a72a9d4186e65bea559
)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(eigen)

# ── CLI11 (header-only argument parsing) ──────────────────────────
FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG        v2.4.1
)
FetchContent_MakeAvailable(cli11)

# ── Pipeline executable ───────────────────────────────────────────
add_executable(clnj_pipeline
    src/main.cpp
    src/fasta_parser.cpp
    src/one_hot.cpp
    src/distance.cpp
    src/kmer_distance.cpp
    src/alignment_utils.cpp
    src/iqtree_interface.cpp
    src/algorithm2.cpp
    src/mlvmst.cpp
    src/distance_oracle.cpp
    src/nj.cpp
    src/mfnj.cpp
    src/bionj.cpp
    src/clnj.cpp
    src/tree_analysis.cpp
    src/mdl_clustering.cpp
    src/tree_state.cpp
    src/online_insertion.cpp
    src/pipeline.cpp
)

target_include_directories(clnj_pipeline PRIVATE include)
target_link_libraries(clnj_pipeline PRIVATE Eigen3::Eigen CLI11::CLI11)

# Release flags for performance
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()
target_compile_options(clnj_pipeline PRIVATE
    $<$<CONFIG:Release>:-O2 -DNDEBUG>
)
