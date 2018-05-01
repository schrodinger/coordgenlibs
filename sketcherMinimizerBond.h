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
#include "CoordgenConfig.hpp"


class sketcherMinimizerRing;
class sketcherMinimizerAtom;


struct sketcherMinimizerBondStereoInfo {
    enum sketcherMinimizerBondStereo {
        cis,
        trans,
        unspecified
    };

    sketcherMinimizerAtom* atom1 = nullptr;
    sketcherMinimizerAtom* atom2 = nullptr;
    sketcherMinimizerBondStereo stereo = unspecified;
};

/*class to represent a covalent bond*/
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
    sketcherMinimizerBond(sketcherMinimizerAtom* at1, sketcherMinimizerAtom*at2) :
    sketcherMinimizerBond() {
        startAtom = at1;
        endAtom = at2;
    }

    virtual ~sketcherMinimizerBond(){};

    virtual bool isResidueInteraction() { return false; }
    sketcherMinimizerAtom* startAtom;
    sketcherMinimizerAtom* endAtom;
    sketcherMinimizerAtom* getStartAtom() const { return startAtom; }
    sketcherMinimizerAtom* getEndAtom() const { return endAtom; }

    /*return bond order*/
    int getBondOrder() const { return bondOrder; }

    void setBondOrder(int order) {bondOrder = order;}


    void setStereoChemistry(sketcherMinimizerBondStereoInfo stereo) {m_stereo = stereo;}

    void EXPORT_COORDGEN setAbsoluteStereoFromStereoInfo();

    /*return true if the bond is part of a small ring (i.e. 8 members or less)*/
    bool isInSmallRing() const;

    /*return true if the bond is part of a macrocycle*/
    bool isInMacrocycle() const;

    /*return true if the bond is to a terminal atom*/
    bool isTerminal() const;

    /*
     does this bond separate two rigid fragments?
     i.e. is bond a single rotatable bond to a non-terminal atom?
     */
    bool EXPORT_COORDGEN isInterFragment() const;
    bool isStereo() const;

    /*given atom1 and atom2 as substituents on the two sides of a double bond,
     should they be put in cis?
     */
    bool markedAsCis(sketcherMinimizerAtom* atom1,
                     sketcherMinimizerAtom* atom2) const;

    /*flip the current bond, mirroring the coordinates of all the atoms on one side of it*/
    void EXPORT_COORDGEN flip();

    /*get the atom bound to the start atom of the bond with the highest CIP priority*/
    sketcherMinimizerAtom* startAtomCIPFirstNeighbor() const;

    /*get the atom bound to the end atom of the bond with the highest CIP priority*/
    sketcherMinimizerAtom* endAtomCIPFirstNeighbor() const;

    /*return true if the E/Z stereochemistry as read from the atoms coordinates matches the label*/
    bool checkStereoChemistry() const;
    int bondOrder;
    bool skip;
    bool isZEActive; // does it have  a Z and E form?
    bool isZ; // used for double bonds to distinguish Z from E form. bonds
              // default to E

    int m_chmN = -1; // idx of the corresponding ChmAtom if molecule comes from 3d


    sketcherMinimizerBondStereoInfo m_stereo;

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
