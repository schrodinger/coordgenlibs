set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.2)
project(coordgen)

file(GLOB SOURCES "*.cpp")
file(GLOB MAEPARSER_SOURCES "maeparser/*.cpp")
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})
add_library(coordgenlibs SHARED ${SOURCES})
add_library(maeparser SHARED ${MAEPARSER_SOURCES})
add_executable(example example/example.cpp)
target_link_libraries(coordgenlibs PUBLIC maeparser)
target_link_libraries(example coordgenlibs)
target_link_libraries(example maeparser)
target_compile_definitions(coordgenlibs PRIVATE IN_COORDGEN)
set_property(TARGET coordgenlibs PROPERTY CXX_VISIBILITY_PRESET "hidden") # used for linux / windows to mark visibility
enable_testing()

