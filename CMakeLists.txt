cmake_minimum_required(VERSION 3.2)
project(coordgen)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Options
option(COORDGEN_RIGOROUS_BUILD "Abort the build if the compiler issues \
       any warnings" OFF )
option(COORDGEN_BUILD_TESTS "Whether test executables should be built" ON)
option(COORDGEN_BUILD_EXAMPLE "Whether to build the sample executable" ON)
option(COORDGEN_BUILD_SHARED_LIBS "Build coordgen as a shared library \
       (turn off for a static one)" ON)

# Use the maeparser_DIR variable to tell CMake where to search for the
# maeparser library.

set(MAEPARSER_VERSION "master" CACHE STRING "maeparser tag to build if \
    a compiled library is not found")

if(MSVC)
    # C4251 disables warnings for export STL containers as arguments
    # (returning a vector of things)
    add_definitions(/wd4251 /wd4275 /wd4996 /D_SCL_SECURE_NO_WARNINGS
                    /D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-DBOOST_ALL_NO_LIB)
endif(MSVC)

if(COORDGEN_RIGOROUS_BUILD)
if(MSVC)
    add_definitions(/WX)
else(MSVC)
    add_definitions(-Wall -Wextra -Werror)
endif(MSVC)
endif(COORDGEN_RIGOROUS_BUILD)


# Source files & headers
file(GLOB SOURCES "*.cpp")

# Build Targets & Configuration -- coordgen library
if(COORDGEN_BUILD_SHARED_LIBS)
    add_library(coordgen SHARED ${SOURCES})
    target_compile_definitions(coordgen PRIVATE "IN_COORDGEN")
    set_property(TARGET coordgen PROPERTY CXX_VISIBILITY_PRESET "hidden")
else(COORDGEN_BUILD_SHARED_LIBS)
    add_library(coordgen STATIC ${SOURCES})
    target_compile_definitions(coordgen PRIVATE "STATIC_MAEPARSER")
    target_compile_definitions(coordgen PRIVATE "STATIC_COORDGEN")
endif(COORDGEN_BUILD_SHARED_LIBS)

# Dependencies
if(TARGET maeparser)
    message(STATUS "Using externally defined maeparser target to "
    "build coordgen")
else()
    include(CoordgenUtils)
    set(MAEPARSER_BUILD_SHARED_LIBS ${COORDGEN_BUILD_SHARED_LIBS} CACHE BOOL "Library type for maeparser")
    find_or_clone_maeparser()
endif()
include_directories(${maeparser_INCLUDE_DIRS})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

target_link_libraries(coordgen ${maeparser_LIBRARIES})

set_target_properties(coordgen
    PROPERTIES
        VERSION 1.3.2
        SOVERSION 1
)

# Install configuration
install(TARGETS coordgen
    EXPORT coordgen-targets
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib)

install(FILES
    CoordgenConfig.hpp
    CoordgenFragmentBuilder.h
    CoordgenFragmenter.h
    CoordgenMacrocycleBuilder.h
    CoordgenMinimizer.h
    sketcherMinimizerAtom.h
    sketcherMinimizerBendInteraction.h
    sketcherMinimizerBond.h
    sketcherMinimizerClashInteraction.h
    sketcherMinimizerEZConstrainInteraction.h
    sketcherMinimizerFragment.h
    sketcherMinimizer.h
    sketcherMinimizerInteraction.h
    sketcherMinimizerMarchingSquares.h
    sketcherMinimizerMaths.h
    sketcherMinimizerMolecule.h
    sketcherMinimizerResidue.h
    sketcherMinimizerResidueInteraction.h
    sketcherMinimizerRing.h
    sketcherMinimizerStretchInteraction.h
    DESTINATION include/coordgen)

install(FILES
    templates.mae
    DESTINATION share/coordgen)

install(EXPORT coordgen-targets
    FILE ${PROJECT_NAME}-config.cmake
    DESTINATION lib/cmake)

# Example
if(COORDGEN_BUILD_EXAMPLE)
    add_subdirectory(example_dir)
endif(COORDGEN_BUILD_EXAMPLE)

# Tests
if(COORDGEN_BUILD_TESTS)
    set(MEMORYCHECK_COMMAND_OPTIONS "--tool=memcheck --time-stamp=yes \
        --num-callers=20 --gen-suppressions=all --leak-check=full \
        --show-reachable=no --trace-children=yes --error-exitcode=29")
    include(CTest)
    add_subdirectory(test)
endif(COORDGEN_BUILD_TESTS)
