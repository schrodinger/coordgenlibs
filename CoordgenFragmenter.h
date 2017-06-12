/*
   Contributors: Nicola Zonta
   Copyright Schrodinger, LLC. All rights reserved
 */

#ifndef COORDGEN_FRAGMENTER_H
#define COORDGEN_FRAGMENTER_H

#include <vector>

class sketcherMinimizerFragment;
class sketcherMinimizerMolecule;
class sketcherMinimizerBond;
class sketcherMinimizerRing;
class sketcherMinimizerAtom;



/*
 class to divide a molecule into rigid fragments
 */

class CoordgenFragmenter
{
  public:
    static void splitIntoFragments(sketcherMinimizerMolecule* molecule);

  private:
    static void
    joinFragments(sketcherMinimizerFragment* fragment1,
                  sketcherMinimizerFragment* fragment2,
                  std::vector<sketcherMinimizerFragment*>& fragments);
    static void processInterFragmentBond(
        sketcherMinimizerBond* bond,
        std::vector<sketcherMinimizerFragment*>& fragments);
    static void processBondInternalToFragment(
        sketcherMinimizerBond* bond,
        std::vector<sketcherMinimizerFragment*>& fragments);
    static void addBondInformation(sketcherMinimizerBond* bond);
    static void addRingInformation(sketcherMinimizerRing* ring);
    static void
    initializeInformation(std::vector<sketcherMinimizerFragment*> fragments,
                          sketcherMinimizerMolecule* molecule);
    static void setChainInfo(sketcherMinimizerFragment* fragment);
    static bool setFixedInfo(sketcherMinimizerFragment* fragment);
    static bool setConstrainedInfo(sketcherMinimizerFragment* fragment);
    static bool isAtomFixed(const sketcherMinimizerAtom* atom);
    static bool isAtomConstrained(const sketcherMinimizerAtom* atom);
    static bool isChain(const sketcherMinimizerFragment* fragment);
    static bool hasPriority(const sketcherMinimizerFragment* fragment1,
                            const sketcherMinimizerFragment* fragment2);
    static int getValueOfCheck(const sketcherMinimizerFragment* fragment,
                               int checkN, bool& checkNoMore);
    static sketcherMinimizerFragment*
    findMainFragment(std::vector<sketcherMinimizerFragment*> fragments);
    static sketcherMinimizerFragment*
    considerChains(std::vector<sketcherMinimizerFragment*> fragments,
                   sketcherMinimizerFragment* mainFragment);
    static unsigned int
    acceptableChainLength(sketcherMinimizerFragment* mainFragment);
    static std::vector<sketcherMinimizerFragment*>
    findLongestChain(std::vector<sketcherMinimizerFragment*> fragments);
    static void addParentRelationsToFragments(
        sketcherMinimizerFragment* mainFragment,
        std::vector<sketcherMinimizerFragment*> fragments);
    static void
    orderFragments(std::vector<sketcherMinimizerFragment*>& fragments,
                   sketcherMinimizerFragment* mainFragment);
};

#endif /* defined(COORDGEN_FRAGMENTER_H) */
