cmake_minimum_required(VERSION 3.20)
project(L1_Kateryna_Ivanchenko)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

find_package(Threads REQUIRED)

# Explicitly tell CMake where Homebrew's libomp is located for AppleClang
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
find_package(OpenMP REQUIRED COMPONENTS CXX)

# Sequential target (only main.cpp)
add_executable(L1_Kateryna_Ivanchenko_Sequential main.cpp)

# Parallel targets
add_executable(L1_Kateryna_Ivanchenko_Parallel L1a.cpp)
add_executable(L1_Kateryna_Ivanchenko_Parallel2 L1b.cpp)

# Copy data files to build folder
configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/data.json"
        "${CMAKE_CURRENT_BINARY_DIR}/data.json"
        COPYONLY
)

configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/data1.json"
        "${CMAKE_CURRENT_BINARY_DIR}/data1.json"
        COPYONLY
)

# Link targets
target_link_libraries(L1_Kateryna_Ivanchenko_Sequential PRIVATE Threads::Threads)
target_link_libraries(L1_Kateryna_Ivanchenko_Parallel PRIVATE Threads::Threads)
target_link_libraries(L1_Kateryna_Ivanchenko_Parallel2 PRIVATE Threads::Threads OpenMP::OpenMP_CXX)