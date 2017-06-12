/*
 *  sketcherMinimizerMolecule.h
 *
 *  Created by Nicola Zonta on 03/05/2010.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */
#include <vector>
#include <iostream>

#ifndef sketcherMINIMIZERMOLECULE_H
#define sketcherMINIMIZERMOLECULE_H


class sketcherMinimizerAtom;
class sketcherMinimizerBond;
class sketcherMinimizerRing;
class sketcherMinimizerPointF;


class sketcherMinimizerFragment;

class  sketcherMinimizerMolecule
{
  public:
    sketcherMinimizerMolecule();
    ~sketcherMinimizerMolecule();


//    void fromChmMol(ChmMol& mol);

    std::vector<sketcherMinimizerAtom*>& getAtoms() { return _atoms; }
    std::vector<sketcherMinimizerBond*>& getBonds() { return _bonds; }
    std::vector<sketcherMinimizerRing*>& getRings() { return _rings; }
    std::vector<sketcherMinimizerFragment*>& getFragments()
    {
        return _fragments;
    }
    void setFragments(std::vector<sketcherMinimizerFragment*> fragments)
    {
        _fragments = fragments;
    }

    void requireMinimization();
    bool minimizationIsRequired();
    std::vector<sketcherMinimizerAtom*> _atoms;
    std::vector<sketcherMinimizerBond*> _bonds;
    std::vector<sketcherMinimizerRing*> _rings;
    std::vector<sketcherMinimizerBond*> m_proximityRelations;

    std::vector<sketcherMinimizerFragment*> _fragments;

    void setMainFragment(sketcherMinimizerFragment* fragment)
    {
        m_mainFragment = fragment;
    }
    sketcherMinimizerFragment* getMainFragment() { return m_mainFragment; }

    bool fixed;
    bool hasFixedFragments, hasConstrainedFragments;
    bool needToAlignNonRingAtoms;
    bool needToAlignWholeMolecule;
    bool isPlaced; // used by arrangeMultipleMolecules

    int totalCharge();
    void boundingBox(sketcherMinimizerPointF& min,
                     sketcherMinimizerPointF& max);
    sketcherMinimizerPointF center();
    static void forceUpdateStruct(std::vector<sketcherMinimizerAtom*>& atoms,
                                  std::vector<sketcherMinimizerBond*>& bonds,
                                  std::vector<sketcherMinimizerRing*>& rings);
    static void
    assignBondsAndNeighbors(std::vector<sketcherMinimizerAtom*>& atoms,
                            std::vector<sketcherMinimizerBond*>& bonds);

    static void findRings(std::vector<sketcherMinimizerBond*>& bonds,
                          std::vector<sketcherMinimizerRing*>& rings);
    static sketcherMinimizerRing* closeRing(sketcherMinimizerBond* bond);
    static void addRing(sketcherMinimizerRing* ring,
                        std::vector<sketcherMinimizerRing*>& rings);

  private:
    sketcherMinimizerFragment* m_mainFragment;
    bool m_requireMinimization;
};

#endif // sketcherMINIMIZERMOLECULE_H
