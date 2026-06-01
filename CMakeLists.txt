cmake_minimum_required(VERSION 3.16)
project(prot2exon VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization flags. NOTE: do NOT bake in -march=native by default — the binary
# ships in pre-built PyPI wheels (cibuildwheel), so it must run on any CPU of the
# target architecture, not just the machine that compiled it. Opt in to
# host-specific tuning with -DPROT2EXON_NATIVE=ON for a local source build.
option(PROT2EXON_NATIVE "Tune the build for the host CPU (-march=native)" OFF)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
if(PROT2EXON_NATIVE)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Find required packages
find_package(Threads REQUIRED)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
set(SOURCES
    src/gtf_parser.cpp
    src/domain_mapper.cpp
    src/output_writer.cpp
    src/utils.cpp
)

# Create executable
add_executable(prot2exon src/main.cpp ${SOURCES})

# Link libraries
target_link_libraries(prot2exon Threads::Threads)

# Optional: Enable OpenMP for parallel processing
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(prot2exon OpenMP::OpenMP_CXX)
    target_compile_definitions(prot2exon PRIVATE USE_OPENMP)
endif()

# When building the PyPI wheel on Linux, statically link the C++ / GCC runtime so
# the bundled executable doesn't depend on the build host's libstdc++ ABI — it
# must run on whatever (older) system the user pip-installs onto.
if(DEFINED SKBUILD AND CMAKE_SYSTEM_NAME STREQUAL "Linux"
   AND CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_link_options(prot2exon PRIVATE -static-libstdc++ -static-libgcc)
endif()

# Install
if(DEFINED SKBUILD)
    # Building a Python wheel (scikit-build-core): bundle the binary INSIDE the
    # package as prot2exon/_bin/prot2exon-core, so a `pip install` ships a ready
    # executable and the `prot2exon` console-script finds it with nothing on PATH.
    install(PROGRAMS $<TARGET_FILE:prot2exon>
            DESTINATION ${SKBUILD_PLATLIB_DIR}/prot2exon/_bin
            RENAME prot2exon-core)
else()
    # Source / conda install: drop it on PATH the usual way.
    install(TARGETS prot2exon DESTINATION bin)
endif()
