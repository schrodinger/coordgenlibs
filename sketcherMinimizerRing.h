/*
 *  sketcherMinimizerRing.h
 *
 *  Created by Nicola Zonta on 03/05/2010.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#include <vector>
#include <iostream>
#include "sketcherMinimizerMaths.h"
#include "CoordgenConfig.hpp"


#ifndef sketcherMINIMIZERRING_H
#define sketcherMINIMIZERRING_H


class sketcherMinimizerAtom;
class sketcherMinimizerPointF;
class sketcherMinimizerBond;


/*class to represent a ring*/
class  sketcherMinimizerRing
{
  public:
    EXPORT_COORDGEN sketcherMinimizerRing();
    EXPORT_COORDGEN ~sketcherMinimizerRing();

    /*rings that share atoms with this*/
    std::vector<sketcherMinimizerRing*> fusedWith;

    /*list of atoms that are shared with each other ring in fusedWith*/
    std::vector<std::vector<sketcherMinimizerAtom*>> fusionAtoms;

    /*list of bonds in the case of two rings attached by a double bond*/
    std::vector<sketcherMinimizerBond*> fusionBonds;
    bool visited, coordinatesGenerated, side /*not central */;
    std::vector<sketcherMinimizerAtom*> getAtoms() const { return _atoms; }
    int size() const { return (int) _atoms.size(); }
    bool isMacrocycle() const { return size() >= MACROCYCLE; }
    std::vector<sketcherMinimizerAtom*> _atoms;
    std::vector<sketcherMinimizerBond*> _bonds;

    /*return the coordinates of the center of the ring*/
    sketcherMinimizerPointF findCenter();

    /*return true if the ring is benzene*/
    bool isBenzene();

    /*return true if the given point is inside the ring*/
    bool contains(sketcherMinimizerPointF p);

    /*return true if the given atom is part of the ring*/
    bool EXPORT_COORDGEN containsAtom(const sketcherMinimizerAtom* a) const;

    /*return true if the given bond is part of the ring*/
    bool containsBond(sketcherMinimizerBond* b);

    /*return true if this is fused with ring*/
    bool isFusedWith(sketcherMinimizerRing* ring);
    std::vector<sketcherMinimizerAtom*>

    /*return the common atoms between this and ring*/
    getFusionAtomsWith(const sketcherMinimizerRing* ring) const;

    /*convenience function used by the SSSR algorithm*/
    bool sameAs(sketcherMinimizerRing* ring);
    bool isAromatic(); // not chemically accurate, but good enough for minimizer
};

#endif // sketcherMINIMIZERRING_H
