set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.2)
project(coordgen)

file(GLOB SOURCES "*.cpp")
file(GLOB MAEPARSER_SOURCES "maeparser/*.cpp")
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})
add_library(maeparser ${MAEPARSER_SOURCES})
add_library(coordgenlibs ${SOURCES})
add_executable(example example/example.cpp)
target_link_libraries(example coordgenlibs)
enable_testing()

