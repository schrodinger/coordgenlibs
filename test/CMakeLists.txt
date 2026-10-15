# Use the example executable as a test
add_test(example ${CMAKE_BINARY_DIR}/example_dir/example)

add_executable(test_coordgen test_coordgen.cpp)

if(COORDGEN_BUILD_SHARED_LIBS)
    target_compile_definitions(test_coordgen PRIVATE "BOOST_ALL_DYN_LINK")
else(COORDGEN_BUILD_SHARED_LIBS)
    set(Boost_USE_STATIC_LIBS ON)
    target_compile_definitions(test_coordgen PRIVATE "STATIC_MAEPARSER")
    target_compile_definitions(test_coordgen PRIVATE "STATIC_COORDGEN")
endif(COORDGEN_BUILD_SHARED_LIBS)

find_package(Boost COMPONENTS filesystem unit_test_framework REQUIRED)

target_link_libraries(test_coordgen coordgen ${maeparser_LIBRARIES} ${boost_link_options} Boost::unit_test_framework Boost::filesystem)

# Set the path for the input files
get_filename_component(TEST_SAMPLES_PATH ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)
target_compile_definitions(test_coordgen PRIVATE "TEST_SAMPLES_PATH=\"${TEST_SAMPLES_PATH}\"")
target_compile_definitions(test_coordgen PRIVATE "SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\"")

add_test(NAME test_coordgen COMMAND ${CMAKE_BINARY_DIR}/test/test_coordgen
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test)


add_executable(test_smilesparser test_smilesparser.cpp)

if(COORDGEN_BUILD_SHARED_LIBS)
    target_compile_definitions(test_smilesparser PRIVATE "BOOST_ALL_DYN_LINK")
else(COORDGEN_BUILD_SHARED_LIBS)
    set(Boost_USE_STATIC_LIBS ON)
    target_compile_definitions(test_smilesparser PRIVATE "STATIC_COORDGEN")
endif(COORDGEN_BUILD_SHARED_LIBS)

target_link_libraries(test_smilesparser coordgen ${boost_link_options} Boost::unit_test_framework Boost::filesystem)

add_test(NAME test_smilesparser COMMAND ${CMAKE_BINARY_DIR}/test/test_coordgen
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test)


# Its own executable: the race it checks exists only before the first template
# load in a process.
add_executable(test_threads test_threads.cpp)

if(COORDGEN_BUILD_SHARED_LIBS)
    target_compile_definitions(test_threads PRIVATE "BOOST_ALL_DYN_LINK")
else(COORDGEN_BUILD_SHARED_LIBS)
    set(Boost_USE_STATIC_LIBS ON)
    target_compile_definitions(test_threads PRIVATE "STATIC_COORDGEN")
endif(COORDGEN_BUILD_SHARED_LIBS)

find_package(Threads REQUIRED)
target_link_libraries(test_threads coordgen ${boost_link_options} Boost::unit_test_framework Threads::Threads)

# One run catches the race only some of the time, so run it in many fresh
# processes: 50 runs are needed to catch it reliably on a 2-core machine.
# Turning off glibc's per-thread malloc cache and scribbling over freed memory
# makes a use of the freed template list crash instead of reading stale
# pointers; without them a 2-core machine almost never fails.
foreach(run RANGE 1 50)
    add_test(NAME test_threads_${run} COMMAND ${CMAKE_BINARY_DIR}/test/test_threads
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test)
    set_tests_properties(test_threads_${run} PROPERTIES ENVIRONMENT
        "MALLOC_PERTURB_=165;GLIBC_TUNABLES=glibc.malloc.tcache_count=0")
endforeach()