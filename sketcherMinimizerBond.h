/*
 *  sketcherMinimizerBond.h
 *
 *  Created by Nicola Zonta on 03/05/2010.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#ifndef sketcherMINIMIZERBOND_H
#define sketcherMINIMIZERBOND_H


#include <cstddef>
#include <vector>

class sketcherMinimizerRing;
class sketcherMinimizerAtom;

class  sketcherMinimizerBond
{
  public:
    sketcherMinimizerBond()
        : startAtom(NULL), endAtom(NULL), bondOrder(1), skip(false),
          isZEActive(false), isZ(false), isWedge(false), isReversed(false),
          hasStereochemistryDisplay(false), _SSSRVisited(false),
          _SSSRParentAtStart(true), m_ignoreZE(false), _SSSRParent(NULL),
          rings()
    {
    }
    virtual ~sketcherMinimizerBond(){};

    virtual bool isResidueInteraction() { return false; }
    sketcherMinimizerAtom* startAtom;
    sketcherMinimizerAtom* endAtom;
    sketcherMinimizerAtom* getStartAtom() const { return startAtom; }
    sketcherMinimizerAtom* getEndAtom() const { return endAtom; }
    int getBondOrder() const { return bondOrder; }
    bool isInSmallRing() const;
    bool isInMacrocycle() const;
    bool isTerminal() const;

    /*
     does this bond separate two rigid fragments?
     i.e. is bond a single rotatable bond to a non-terminal atom?
     */
    bool isInterFragment() const;
    bool isStereo() const;
    /*given atom1 and atom2 as substituents on the two sides of a double bond,
     should they be put in cis?
     */
    bool markedAsCis(sketcherMinimizerAtom* atom1,
                     sketcherMinimizerAtom* atom2) const;

    void flip();
    sketcherMinimizerAtom* startAtomCIPFirstNeighbor() const;
    sketcherMinimizerAtom* endAtomCIPFirstNeighbor() const;
    bool checkStereoChemistry() const;
    int bondOrder;
    bool skip;
    bool isZEActive; // does it have  a Z and E form?
    bool isZ; // used for double bonds to distinguish Z from E form. bonds
              // default to E
    bool isWedge;
    bool isReversed;
    bool hasStereochemistryDisplay;

    bool _SSSRVisited;
    bool _SSSRParentAtStart;
    bool m_ignoreZE;
    sketcherMinimizerBond* _SSSRParent;
    std::vector<sketcherMinimizerRing*> rings;
};

#endif // sketcherMINIMIZERBOND_H
