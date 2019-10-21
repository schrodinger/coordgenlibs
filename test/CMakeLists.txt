# Use the example executable as a test
add_test(example ${CMAKE_BINARY_DIR}/example_dir/example)

# Require Boost testing framework for real tests
find_package(Boost COMPONENTS unit_test_framework filesystem REQUIRED)

add_executable(test_coordgen test_coordgen.cpp)
target_link_libraries(test_coordgen coordgen Boost::unit_test_framework Boost::filesystem)

# Set the path for the input files
get_filename_component(TEST_SAMPLES_PATH ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)
target_compile_definitions(test_coordgen PRIVATE "TEST_SAMPLES_PATH=\"${TEST_SAMPLES_PATH}\"")

add_test(NAME test_coordgen COMMAND ${CMAKE_BINARY_DIR}/test/test_coordgen
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test)
