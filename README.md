# CoordgenLibs

[![Build Status](https://travis-ci.org/schrodinger/coordgenlibs.svg?branch=master)](https://travis-ci.org/schrodinger/coordgenlibs)
[![Build_Status](https://ci.appveyor.com/api/projects/status/github/schrodinger/coordgenlibs?branch=master&svg=true)](https://ci.appveyor.com/project/torcolvin/coordgenlibs-3h7cs)


This is **Schrodinger, Inc's** 2D coordinate generation.  It was formerly proprietary code, but is now released under the [BSD license](https://github.com/schrodinger/coordgenlibs/blob/master/LICENSE).  The emphasis of these algorithms are on quality of 2D coordinates rather than speed of generation.  The algorithm distinguishes itself from many others by doing well with both macrocycles and metal complexes.  It also does extremely well on typical drug-like small molecules, and has been validated on millions of compounds.

Schrodinger intends to continue to contribute to this code as it still uses it inside its products, but will also be happy if others contribute pull-requests when there are improvements they would like to make.  We'll also be happy to hear bug reports or feature requests from use of this code, though make no guarantee on our ability to process these.

### Documentation
Examples and documentation will be added/improved over time

### Usage example
Code for a sample executable is provided in the `example_dir` directory. Building the example executable is enabled by default, but can be disabled by means of the `COORDGEN_BUILD_EXAMPLE` option.

### Automated Testing
Automated testing is still primarily taking place inside Schrodinger's internal build system, although tests are incrementally being added to the `testing` directory. Building the tests is enabled by default, but can be disabled by means of the `COORDGEN_BUILD_TESTS` option.

Memory debugging is, by default, configured to use `valgrind`. It can be run on the tests by passing `-DCMAKE_BUILD_TYPE=Debug` to cmake, to enable building the debugging symbols, and then using `ctest -T memcheck` inside the build directory.

