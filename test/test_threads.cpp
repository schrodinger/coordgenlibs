#define BOOST_TEST_MODULE Test_Coordgen_Threads

#include <boost/test/unit_test.hpp>

#include <atomic>
#include <thread>
#include <vector>

#include "../sketcherMinimizer.h"
#include "coordgenBasicSMILES.h"

using namespace schrodinger;

// sketcherMinimizer::m_templates is loaded on first use and shared by every
// thread in the process. This is its own executable so that the first use is
// here: once anything has loaded the templates, there is nothing left to race.
BOOST_AUTO_TEST_CASE(ConcurrentFirstTemplateLoad)
{
    BOOST_REQUIRE(sketcherMinimizer::m_templates.getTemplates().empty());

    // A 21-membered ring closed through the para positions of a benzene. Two
    // fused rings that cannot be reduced to one central ring, so layout
    // searches the templates.
    const std::string smiles = "C1CCCCCC2=CC=C(C=C2)CCCCCCCCCCC1";
    const int nThreads = 16;

    std::atomic<int> ready(0);
    std::atomic<bool> go(false);
    std::vector<std::thread> threads;
    for (int i = 0; i < nThreads; ++i) {
        threads.emplace_back([&] {
            auto* mol = approxSmilesParse(smiles);
            sketcherMinimizer minimizer;
            ++ready;
            while (!go) {
                std::this_thread::yield();
            }
            minimizer.initialize(mol); // minimizer takes ownership of mol
            minimizer.runGenerateCoordinates();
        });
    }
    while (ready < nThreads) {
        std::this_thread::yield();
    }
    go = true;
    for (auto& t : threads) {
        t.join();
    }

    // Without a guarded load this test usually crashes before it gets here,
    // and always reports under ThreadSanitizer or AddressSanitizer.
    BOOST_CHECK(!sketcherMinimizer::m_templates.getTemplates().empty());
}
