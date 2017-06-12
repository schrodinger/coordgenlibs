/*
 *
 *
 *  Created by Nicola Zonta on 03/05/2010.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#ifndef sketcherMINIMIZERFRAGMENT
#define sketcherMINIMIZERFRAGMENT

#include <vector>
#include <map>
#include <assert.h>
#include <iostream>
#include "sketcherMinimizerMaths.h"

class sketcherMinimizerAtom;
class sketcherMinimizerBond;
class sketcherMinimizerRing;
class sketcherMinimizerFragment;

/*
 abstract class for fragment degree of freedom
 */
class CoordgenFragmentDOF
{
  public:
    CoordgenFragmentDOF(sketcherMinimizerFragment* fragment);
    virtual ~CoordgenFragmentDOF();
    void storeCurrentValueAsOptimal();
    void setToOptimalValue();
    void changeState();
    void addAtom(sketcherMinimizerAtom* atom);
    virtual float getCurrentPenalty() const;
    virtual int numberOfStates() const = 0;
    virtual int tier() const = 0;
    virtual void apply() const = 0;
    sketcherMinimizerFragment* getFragment() const;
    short unsigned int getCurrentState();
    void setState(short unsigned int state);
    short unsigned int m_currentState, m_optimalState;

  protected:
    std::vector<sketcherMinimizerAtom*> m_atoms;
    sketcherMinimizerFragment* m_fragment;
};

class CoordgenRotateFragmentDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenRotateFragmentDOF(sketcherMinimizerFragment* fragment);
    int numberOfStates() const;
    int tier() const;
    void apply() const;
    float getCurrentPenalty() const;
};

class CoordgenFlipFragmentDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenFlipFragmentDOF(sketcherMinimizerFragment* fragment);
    int numberOfStates() const;
    int tier() const;
    void apply() const;
    float getCurrentPenalty() const;
};

class CoordgenScaleAtomsDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenScaleAtomsDOF(sketcherMinimizerAtom* pivotAtom);
    int numberOfStates() const;
    int tier() const;
    void apply() const;
    float getCurrentPenalty() const;

  private:
    sketcherMinimizerAtom* m_pivotAtom;
};

class CoordgenScaleFragmentDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenScaleFragmentDOF(sketcherMinimizerFragment* fragment);
    int numberOfStates() const;
    int tier() const;
    void apply() const;
    float getCurrentPenalty() const;
};

class CoordgenChangeParentBondLengthFragmentDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenChangeParentBondLengthFragmentDOF(
        sketcherMinimizerFragment* fragment);
    int numberOfStates() const;
    int tier() const;
    void apply() const;
    float getCurrentPenalty() const;
};

class CoordgenInvertBondDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenInvertBondDOF(sketcherMinimizerAtom* pivotAtom,
                          sketcherMinimizerAtom* boundAtom);
    int numberOfStates() const;
    int tier() const;
    void apply() const;
    float getCurrentPenalty() const;

  private:
    sketcherMinimizerAtom* m_pivotAtom;
    sketcherMinimizerAtom* m_boundAtom;
};

class CoordgenFlipRingDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenFlipRingDOF(sketcherMinimizerRing* ring,
                        std::vector<sketcherMinimizerAtom*> fusionAtoms);
    int numberOfStates() const;
    int tier() const;
    void apply() const;
    float getCurrentPenalty() const;

  private:
    sketcherMinimizerAtom* m_pivotAtom1;
    sketcherMinimizerAtom* m_pivotAtom2;
    int m_penalty;
};

class sketcherMinimizerFragment
{
  public:
    sketcherMinimizerFragment();
    ~sketcherMinimizerFragment();
    unsigned int totalWeight() const;
    unsigned int countDoubleBonds() const;
    unsigned int countHeavyAtoms() const;

    void addAtom(sketcherMinimizerAtom* atom);
    void addBond(sketcherMinimizerBond* bond);
    void addRing(sketcherMinimizerRing* ring);
    void addDof(CoordgenFragmentDOF* dof);
    std::vector<CoordgenFragmentDOF*> getDofs();

    void addInterFragmentBond(sketcherMinimizerBond* bond)
    {
        _interFragmentBonds.push_back(bond);
    }

    std::vector<sketcherMinimizerAtom*> getAtoms() const { return m_atoms; }
    std::vector<sketcherMinimizerBond*> getBonds() const { return m_bonds; }
    std::vector<sketcherMinimizerRing*> getRings() const { return m_rings; }

    std::vector<sketcherMinimizerAtom*>& atoms() { return m_atoms; }
    std::vector<sketcherMinimizerBond*>& bonds() { return m_bonds; }
    std::vector<sketcherMinimizerRing*>& rings() { return m_rings; }
    sketcherMinimizerFragment* getParent() const { return m_parent; }
    void setParent(sketcherMinimizerFragment* parent) { m_parent = parent; }
    void addDofToAtom(sketcherMinimizerAtom* atom, CoordgenFragmentDOF* dof)
    {
        m_dofsForAtom[atom].push_back(dof);
    }
    std::vector<CoordgenFragmentDOF*>&
    getDofsOfAtom(sketcherMinimizerAtom* atom)
    {
        return m_dofsForAtom[atom];
    }
    void setAllCoordinatesToTemplate();
    void storeCoordinateInformation();

    std::vector<sketcherMinimizerBond*> _interFragmentBonds;
    std::vector<sketcherMinimizerFragment*> _children;
    std::map<sketcherMinimizerAtom*, sketcherMinimizerPointF> _coordinates;
    sketcherMinimizerPointF _bondToParentCoordinatesStart;
    sketcherMinimizerPointF _bondToParentCoordinatesEnd;
    bool fixed, isTemplated, constrained;
    bool isChain;
    sketcherMinimizerBond* _bondToParent;
    float longestChainFromHere;
    int numberOfChildrenAtoms;
    float numberOfChildrenAtomsRank;

    void setCoordinates(sketcherMinimizerPointF position, float angle);
    CoordgenFragmentDOF* getFlipDof() const { return m_dofs[0]; }
  private:
    sketcherMinimizerFragment* m_parent;
    std::vector<sketcherMinimizerAtom*> m_atoms;
    std::vector<sketcherMinimizerBond*> m_bonds;
    std::vector<sketcherMinimizerRing*> m_rings;
    std::vector<CoordgenFragmentDOF*> m_dofs;
    std::map<sketcherMinimizerAtom*, std::vector<CoordgenFragmentDOF*>>
        m_dofsForAtom;
};

#endif // sketcherMINIMIZERFRAGMENT
