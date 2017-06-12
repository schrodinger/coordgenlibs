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

#ifndef sketcherMINIMIZERRING_H
#define sketcherMINIMIZERRING_H


class sketcherMinimizerAtom;
class sketcherMinimizerPointF;
class sketcherMinimizerBond;

class  sketcherMinimizerRing
{
  public:
    sketcherMinimizerRing();
    ~sketcherMinimizerRing();
    std::vector<sketcherMinimizerRing*> fusedWith;
    std::vector<std::vector<sketcherMinimizerAtom*>> fusionAtoms;
    std::vector<sketcherMinimizerBond*> fusionBonds;
    bool visited, coordinatesGenerated, side /*not central */;
    std::vector<sketcherMinimizerAtom*> getAtoms() const { return _atoms; }
    int size() const { return (int) _atoms.size(); }
    bool isMacrocycle() const { return size() >= MACROCYCLE; }
    std::vector<sketcherMinimizerAtom*> _atoms;
    std::vector<sketcherMinimizerBond*> _bonds;

    sketcherMinimizerPointF findCenter();
    bool isBenzene();
    bool contains(sketcherMinimizerPointF p);
    bool containsAtom(const sketcherMinimizerAtom* a) const;
    bool containsBond(sketcherMinimizerBond* b);
    bool isFusedWith(sketcherMinimizerRing* ring);
    std::vector<sketcherMinimizerAtom*>
    getFusionAtomsWith(const sketcherMinimizerRing* ring) const;
    bool sameAs(sketcherMinimizerRing* ring);
    bool isAromatic(); // not chemically accurate, but good enough for minimizer
};

#endif // sketcherMINIMIZERRING_H
