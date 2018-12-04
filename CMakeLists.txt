set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.2)
project(coordgen)


file(GLOB SOURCES "*.cpp")
find_package(Boost COMPONENTS iostreams REQUIRED)
find_package(maeparser)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${maeparser_INCLUDE_DIRS})
add_library(coordgen SHARED ${SOURCES})
target_compile_definitions(coordgen PRIVATE IN_COORDGEN)
set_property(TARGET coordgen PROPERTY CXX_VISIBILITY_PRESET "hidden")
SET_TARGET_PROPERTIES(coordgen
    PROPERTIES
        VERSION 1.2
        SOVERSION 1
)
add_executable(example example_dir/example.cpp)

target_link_libraries(coordgen maeparser)
target_link_libraries(example coordgen)
enable_testing()

install(TARGETS coordgen
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
