cmake_minimum_required(VERSION 3.16)

project(wcns_v2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- MPI ----
find_package(MPI REQUIRED)

# ---- CGNS library (local installation) ----
set(CGNS_ROOT "/home/huangwz_lab/library/cgns")
set(CGNS_INCLUDE_DIR "${CGNS_ROOT}/include")
set(CGNS_LIBRARY     "${CGNS_ROOT}/lib/libcgns.so")

if(NOT EXISTS "${CGNS_INCLUDE_DIR}/cgnslib.h")
    message(FATAL_ERROR "CGNS header not found at ${CGNS_INCLUDE_DIR}/cgnslib.h")
endif()

# ---- Source files ----
set(SOURCES
    src/main.cpp
    src/io/cgns_reader.cpp
    src/io/config_reader.cpp
)

add_executable(wcns_v2 ${SOURCES})

target_include_directories(wcns_v2 PRIVATE
    include
    ${CGNS_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_DIRS}
)

target_link_libraries(wcns_v2 PRIVATE ${CGNS_LIBRARY} ${MPI_CXX_LIBRARIES})

# ---- Compiler flags ----
target_compile_options(wcns_v2 PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
)

if(CMAKE_BUILD_TYPE STREQUAL "Release" OR NOT CMAKE_BUILD_TYPE)
    target_compile_options(wcns_v2 PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-O3 -march=native>
    )
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(wcns_v2 PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-g -O0>
    )
endif()
