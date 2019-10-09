# CoordgenLibs

[![Build Status](https://travis-ci.org/schrodinger/coordgenlibs.svg?branch=master)](https://travis-ci.org/schrodinger/coordgenlibs)
[![Build_Status](https://ci.appveyor.com/api/projects/status/github/schrodinger/coordgenlibs?branch=master&svg=true)](https://ci.appveyor.com/project/torcolvin/coordgenlibs-3h7cs)

This is **Schrodinger, Inc's** 2D coordinate generation.  It was formerly proprietary code, but is now released under the [BSD license](https://github.com/schrodinger/coordgenlibs/blob/master/LICENSE).  The emphasis of these algorithms are on quality of 2D coordinates rather than speed of generation.  The algorithm distinguishes itself from many others by doing well with both macrocycles and metal complexes.  It also does extremely well on typical drug-like small molecules, and has been validated on millions of compounds.

Schrodinger intends to continue to contribute to this code as it still uses it inside its products, but will also be happy if others contribute pull-requests when there are improvements they would like to make.  We'll also be happy to hear bug reports or feature requests from use of this code, though make no guarantee on our ability to process these.

## Documentation

Examples and documentation will be added/improved over time

## Usage example

Code for a sample executable is provided in the `example_dir` directory. Building the example executable is enabled by default, but can be disabled by means of the `COORDGEN_BUILD_EXAMPLE` option.

## Automated Testing

Automated testing is still primarily taking place inside Schrodinger's internal build system, although tests are incrementally being added to the `testing` directory. Building the tests is enabled by default, but can be disabled by means of the `COORDGEN_BUILD_TESTS` option.
d
Memory debugging is, by default, configured to use `valgrind`. It can be run on the tests by passing `-DCMAKE_BUILD_TYPE=Debug` to cmake, to enable building the debugging symbols, and then using `ctest -T memcheck` inside the build directory.

## Building from source

### Requirements

To build coordgen, you will need to have the following installed in your system:

- **CMake** version 3.2 or later.
- The development files for the **Boost libraries**. At least the **iostreams** and **regex** components are required. In case of also building the unit tests, the **filesystems** and **unit_test_framework** components will also be required.
- A **C++ compiler** supporting the C++11 standard.
- A compiled instance of the **maeparser library** or its source code.

In case **maeparser** is not available on your system, neither as a compiled library or as source code,if a working `git` executable and an internet connection are available, the builder can automatically download the source and build **maeparser** for you.

### Building

1. Create a build directory inside the the one that contains Coordgen, and move into it:

```bash
mkdir build
cd build
```

1. Run `cmake`, passing the path to the directory where the sources are located (just `..` if you created `build` inside the sources directory). At this point, you should add any required flags to the `cmake` command. Check the 'Options' section in CMakeLists.txt to see which options are available.

```bash
cmake .. -Dmaeparser_DIR=/home/schrodinger/maeparser_install -DCMAKE_INSTALL_PREFIX=/home/schrodinger/coordgen_install`
```

A few notes on the `maeparser_DIR` option:

- CMake will, by default, search your system's default library paths for the maeparser library. If a `CMAKE_INSTALL_PREFIX` was specified, CMake will also search for maeparser there.

- If you used the `CMAKE_INSTALL_PREFIX` to build and install maeparser, you should give the exact same path to `maeparser_DIR`.

- CMake will look for the required maeparser headers and library under the indicated path. In case your headers and library live under different paths, point `maeparser_DIR` to the path where your library lives, and point `maeparser_INCLUDE_DIRS` to the path above your `maeparser/Reader.hpp` header.

- If `maeparser_DIR` was passed to CMake, and the library was not found, CMake will **NOT** download the sources from GitHub (since we expected to find a compiled library).

- Even if the sources cannot be cloned/updated from GitHub, if a copy of maeparser's source is found in the right place (a `maeparser` directory inside Coordgen's source directory), it will be built if no compiled library is available.

1. Build and install:

```bash
make -j install
```
