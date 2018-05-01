set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.2)
project(coordgen)

#include(ExternalProject)

#ExternalProject_Add(
#    maeparser
#    SOURCE_DIR ${coordgen_SOURCE_DIR}/maeparser
#)


file(GLOB SOURCES "*.cpp")
find_package(Boost)
find_package(maeparser)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${maeparser_INCLUDE_DIRS})
add_library(coordgenlibs SHARED ${SOURCES})
add_executable(example example/example.cpp)

target_link_libraries(coordgenlibs maeparser)
target_link_libraries(example coordgenlibs)
enable_testing()

install(TARGETS coordgenlibs
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib)

install(FILES
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
