#define BOOST_TEST_MODULE Test_Coordgen

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include "../sketcherMinimizer.h"
#include "../sketcherMinimizerMaths.h"

#include "maeparser/MaeConstants.hpp"
#include "maeparser/Reader.hpp"

using namespace schrodinger;

const boost::filesystem::path test_samples_path(TEST_SAMPLES_PATH);

namespace
{
std::map<sketcherMinimizerAtom*, int>
getReportingIndices(sketcherMinimizerMolecule& mol)
{
    std::map<sketcherMinimizerAtom*, int> fakeIndices;
    int index = 0;
    for (auto& atom : mol.getAtoms()) {
        fakeIndices.emplace(atom, ++index);
    }
    return fakeIndices;
}

bool areBondsNearIdeal(sketcherMinimizerMolecule& mol,
                       std::map<sketcherMinimizerAtom*, int>& indices)
{
    const float targetBondLength = BONDLENGTH * BONDLENGTH;
    const auto tolerance = static_cast<float>(targetBondLength * 0.1);

    bool passed = true;
    for (auto& bond : mol.getBonds()) {
        auto& startCoordinates = bond->getStartAtom()->getCoordinates();
        auto& endCoordinates = bond->getEndAtom()->getCoordinates();

        const auto sq_distance = sketcherMinimizerMaths::squaredDistance(
            startCoordinates, endCoordinates);
        const auto deviation = sq_distance - targetBondLength;
        if (deviation < -tolerance || deviation > tolerance) {
            std::cerr << "Bond" << indices[bond->getStartAtom()] << '-'
                      << indices[bond->getEndAtom()] << " has length "
                      << sq_distance << " (" << targetBondLength << ")\n";
            passed = false;
        }
    }

    return passed;
}
}

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

    auto indices = getReportingIndices(*mol);
    BOOST_CHECK(areBondsNearIdeal(*mol, indices));
}
