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

    for (auto& atom : mol->getAtoms()) {
        auto c = atom->getCoordinates();

        // This is best we can do with checking: there are precision and
        // rounding issues depending on platform and environment.
        BOOST_CHECK(c.x() != 0 || c.y() != 0);
    }
}
