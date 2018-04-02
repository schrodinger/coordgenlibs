/*
 *  sketcherMinimizerRing.cpp
 *
 *  Created by Nicola Zonta on 24/05/2011.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#include "sketcherMinimizerRing.h"
#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerMaths.h"

using namespace std;

sketcherMinimizerRing::sketcherMinimizerRing()
    : visited(false), coordinatesGenerated(false)
{

    //    assert (0);
    side = false;
}

sketcherMinimizerRing::~sketcherMinimizerRing()
{
}

sketcherMinimizerPointF sketcherMinimizerRing::findCenter()
{
    sketcherMinimizerPointF o(0.f, 0.f);
    for (unsigned int i = 0; i < _atoms.size(); i++)
        o += _atoms[i]->coordinates;
    o /= _atoms.size();
    return o;
}

bool sketcherMinimizerRing::isBenzene()
{
    if (_atoms.size() != 6)
        return false;
    for (unsigned int i = 0; i < _atoms.size(); i++)
        if (_atoms[i]->atomicNumber != 6)
            return false;
    for (unsigned int i = 0; i < _atoms.size(); i++) {
        sketcherMinimizerAtom* a = _atoms[i];
        bool found = false;
        for (unsigned int j = 0; j < a->bonds.size(); j++) {
            if (a->bonds[j]->bondOrder == 2) {
                found = true;
                break;
            }
        }
        if (!found)
            return false;
    }

    return true;
}

bool sketcherMinimizerRing::isAromatic() // not chemically accurate, but good
                                         // enough for minimizer
{
    int bonds = _bonds.size();
    int doubleBonds = 0;
    int NSOCount = 0;
    for (unsigned int i = 0; i < _bonds.size(); i++) {
        if (_bonds[i]->bondOrder == 2)
            doubleBonds++;
    }
    for (unsigned int i = 0; i < _atoms.size(); i++) {
        int an = _atoms[i]->atomicNumber;
        bool doubleBound = false;
        for (unsigned int b = 0; b < _atoms[i]->bonds.size(); b++)
            if (_atoms[i]->bonds[b]->bondOrder == 2)
                doubleBound = true;

        if (!doubleBound)
            if (an == 8 || an == 7 || an == 16)
                NSOCount++;
    }
    if (bonds == 6 && doubleBonds == 3)
        return true;
    if (bonds == 5 && doubleBonds == 2 && NSOCount == 1)
        return true;
    return false;
}

bool sketcherMinimizerRing::containsAtom(const sketcherMinimizerAtom* a) const
{
    for (unsigned int i = 0; i < _atoms.size(); i++)
        if (_atoms[i] == a)
            return true;
    return false;
}
bool sketcherMinimizerRing::containsBond(sketcherMinimizerBond* b)
{
    for (unsigned int i = 0; i < _bonds.size(); i++)
        if (_bonds[i] == b)
            return true;
    return false;
}
bool sketcherMinimizerRing::isFusedWith(sketcherMinimizerRing* ring)
{
    for (unsigned int i = 0; i < fusedWith.size(); i++) {
        if (fusedWith[i] == ring) {
            return true;
        }
    }
    return false;
}
std::vector<sketcherMinimizerAtom*> sketcherMinimizerRing::getFusionAtomsWith(
    const sketcherMinimizerRing* ring) const
{
    for (unsigned int i = 0; i < fusedWith.size(); i++) {
        if (fusedWith[i] == ring) {
            return fusionAtoms[i];
        }
    }
    std::vector<sketcherMinimizerAtom*> empty;
    return empty;
}

bool sketcherMinimizerRing::sameAs(sketcherMinimizerRing* ring)
{
    if (!(_bonds.size() == ring->_bonds.size()))
        return false;
    for (unsigned int i = 0; i < _bonds.size(); i++) {
        if (!ring->containsBond(_bonds[i])) {
            return false;
        }
    }
    return true;
}

bool sketcherMinimizerRing::contains(sketcherMinimizerPointF p)
{

    int n = 0;
    for (unsigned int i = 0; i < _bonds.size(); i++) {
        sketcherMinimizerBond* b = _bonds[i];
        if ((p.y() < b->startAtom->coordinates.y() &&
             p.y() > b->endAtom->coordinates.y()) ||
            (p.y() > b->startAtom->coordinates.y() &&
             p.y() < b->endAtom->coordinates.y())) {
            sketcherMinimizerPointF v =
                b->endAtom->coordinates - b->startAtom->coordinates;
            if (v.y() > SKETCHER_EPSILON || v.y() < -SKETCHER_EPSILON) {
                v *= (p.y() - b->startAtom->coordinates.y()) / v.y();
                v += b->startAtom->coordinates;
                if (p.x() > v.x())
                    n++;
            }
        }
    }
    if ((n % 2) != 0)
        return true;
    return false;
}

//    std::vector <sketcherMinimizerBond *> _bonds;
