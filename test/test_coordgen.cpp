#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test_Coordgen

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include "../sketcherMinimizer.h"

#include "maeparser/MaeConstants.hpp"
#include "maeparser/Reader.hpp"

using namespace schrodinger;

const boost::filesystem::path test_samples_path(TEST_SAMPLES_PATH);

BOOST_AUTO_TEST_CASE(SampleTest)
{
    // A small sample test showing how to import a molecule from a .mae file and
    // initialize the minimizer with it.

    const std::string testfile = (test_samples_path / "test.mae").string();
    const std::vector<std::pair<double, double>> ref_coords(
        {{6.67, -6.66},     {-36.35, -32.58},  {6.88, 42.66},
         {49.93, -32.72},   {-37.25, -82.67},  {-81.13, -7.44},
         {50.61, 68.04},    {-42.27, 42.42},   {-18.43, 85.36},
         {50.93, -82.83},   {93.65, -7.44},    {6.98, -107.07},
         {-80.09, -107.45}, {-124, -33.2},     {-82.25, 42.75},
         {51.42, 118.71},   {6.98, -157.58},   {-123.42, -82.49},
         {101.37, 117.85},  {52.41, 168.7},    {1.78, 119.97},
         {-43.03, -157.56}, {-18.04, -200.86}, {50.27, -182.58},
         {127.06, 160.75},  {96.16, 192.87}});

    mae::Reader r(testfile);
    auto block = r.next(mae::CT_BLOCK);
    BOOST_REQUIRE(block != nullptr);

    auto* mol = mol_from_mae_block(*block);
    BOOST_REQUIRE(mol != nullptr);
    BOOST_REQUIRE_EQUAL(mol->getAtoms().size(), 26);
    BOOST_REQUIRE_EQUAL(mol->getBonds().size(), 26);

    sketcherMinimizer minimizer;
    minimizer.initialize(mol); // minimizer takes ownership of mol
    minimizer.runGenerateCoordinates();

    auto ref_c = ref_coords.begin();
    for (auto& atom : mol->getAtoms()) {
        auto c = atom->getCoordinates();

        BOOST_CHECK_EQUAL(static_cast<float>(ref_c->first), c.x());
        BOOST_CHECK_EQUAL(static_cast<float>(ref_c->second), c.y());

        ++ref_c;
    }
}
