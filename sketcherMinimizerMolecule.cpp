/*
 *  sketcherMinimizerMolecule.cpp
 *
 *  Created by Nicola Zonta on 24/05/2011.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#include "sketcherMinimizerMolecule.h"
#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerMaths.h"
#include "sketcherMinimizerRing.h"
#include "sketcherMinimizerBond.h"
#include <queue>


using namespace std;

sketcherMinimizerMolecule::sketcherMinimizerMolecule()
    : fixed(false),

      hasFixedFragments(false), hasConstrainedFragments(false),
      needToAlignNonRingAtoms(false), needToAlignWholeMolecule(false),
      isPlaced(false), m_mainFragment(NULL), m_requireMinimization(false){};

sketcherMinimizerMolecule::~sketcherMinimizerMolecule()
{
    for (auto ring : _rings) {
        delete ring;
    }
};


sketcherMinimizerAtom* sketcherMinimizerMolecule::addNewAtom()
{
    auto atom = new sketcherMinimizerAtom();
    _atoms.push_back(atom);
    atom->molecule = this;
    return atom;
}


sketcherMinimizerBond* sketcherMinimizerMolecule::addNewBond(sketcherMinimizerAtom* at1,
                                  sketcherMinimizerAtom* at2)
{
    auto bond = new sketcherMinimizerBond(at1, at2);
    _bonds.push_back(bond);
    return bond;
}



int sketcherMinimizerMolecule::totalCharge()
{
    int charge = 0;
    for (unsigned int i = 0; i < _atoms.size(); i++)
        charge += _atoms[i]->charge;
    return charge;
}

void sketcherMinimizerMolecule::boundingBox(sketcherMinimizerPointF& min,
                                            sketcherMinimizerPointF& max)
{
    min.setX(0.f);
    min.setY(0.f);
    max.setX(0.f);
    max.setY(0.f);
    if (_atoms.size()) {
        min = _atoms[0]->coordinates;
        max = _atoms[0]->coordinates;
        for (unsigned int i = 0; i < _atoms.size(); i++) {
            sketcherMinimizerAtom* a = _atoms[i];
            if (a->coordinates.x() < min.x())
                min.setX(a->coordinates.x());
            if (a->coordinates.y() < min.y())
                min.setY(a->coordinates.y());
            if (a->coordinates.x() > max.x())
                max.setX(a->coordinates.x());
            if (a->coordinates.y() > max.y())
                max.setY(a->coordinates.y());
        }
    }
}

void sketcherMinimizerMolecule::requireMinimization()
{
    m_requireMinimization = true;
}

bool sketcherMinimizerMolecule::minimizationIsRequired()
{
    return m_requireMinimization;
}

sketcherMinimizerPointF sketcherMinimizerMolecule::center()
{
    if (!_atoms.size())
        return sketcherMinimizerPointF(0.f, 0.f);
    sketcherMinimizerPointF c(.0f, .0f);
    for (unsigned int i = 0; i < _atoms.size(); i++) {
        c += _atoms[i]->coordinates;
    }
    return c / _atoms.size();
}

void sketcherMinimizerMolecule::assignBondsAndNeighbors(
    std::vector<sketcherMinimizerAtom*>& atoms,
    std::vector<sketcherMinimizerBond*>& bonds)
{
    for (unsigned int i = 0; i < atoms.size(); i++) {
        sketcherMinimizerAtom* atom = atoms[i];
        atom->bonds.clear();
        atom->neighbors.clear();
        atom->residueInteractionPartners.clear();
        atom->residueInteractions.clear();
        atom->rings.clear();
    }

    for (unsigned int i = 0; i < bonds.size(); i++) {

        sketcherMinimizerBond* bond = bonds[i];
        bond->rings.clear();

        if (!bond->isResidueInteraction()) {
            //    bond->_rings.clear ();
            bond->startAtom->bonds.push_back(bond);

            bond->endAtom->neighbors.push_back(bond->startAtom);

            bond->endAtom->bonds.push_back(bond);

            bond->startAtom->neighbors.push_back(bond->endAtom);
        } else {
            bond->startAtom->residueInteractions.push_back(bond);

            bond->endAtom->residueInteractionPartners.push_back(
                bond->startAtom);

            bond->endAtom->residueInteractions.push_back(bond);

            bond->startAtom->residueInteractionPartners.push_back(
                bond->endAtom);
        }
    }
    for (unsigned int i = 0; i < atoms.size(); i++) {
        if (atoms[i]->_implicitHs == -1) atoms[i]->_implicitHs = atoms[i]->findHsNumber();
    }
}

void sketcherMinimizerMolecule::forceUpdateStruct(
    std::vector<sketcherMinimizerAtom*>& atoms,
    std::vector<sketcherMinimizerBond*>& bonds,
    std::vector<sketcherMinimizerRing*>& rings)
{
    sketcherMinimizerMolecule::assignBondsAndNeighbors(atoms, bonds);

    findRings(bonds, rings);
    for (unsigned int i = 0; i < bonds.size(); i++) {
        for (unsigned int j = 0; j < bonds[i]->rings.size(); j++) {
            sketcherMinimizerRing* ring = bonds[i]->rings[j];
            bool found = false;
            for (unsigned int k = 0; k < bonds[i]->startAtom->rings.size();
                 k++) {
                if (bonds[i]->startAtom->rings[k] == ring) {
                    found = true;
                    break;
                }
            }
            if (!found)
                bonds[i]->startAtom->rings.push_back(ring);
            found = false;
            for (unsigned int k = 0; k < bonds[i]->endAtom->rings.size(); k++) {
                if (bonds[i]->endAtom->rings[k] == ring) {
                    found = true;
                    break;
                }
            }
            if (!found)
                bonds[i]->endAtom->rings.push_back(ring);
        }
    }

    for (unsigned int i = 0; i < atoms.size(); i++) {
        for (unsigned int j = 0; j < atoms[i]->rings.size(); j++) {
            atoms[i]->rings[j]->_atoms.push_back(atoms[i]);
        }
    }
}

void sketcherMinimizerMolecule::findRings(
    std::vector<sketcherMinimizerBond*>& bonds,
    std::vector<sketcherMinimizerRing*>& rings)
{
    for (unsigned int i = 0; i < rings.size(); i++)
        delete rings[i];
    rings.clear();
    for (unsigned int i = 0; i < bonds.size(); i++) {
        for (unsigned int ii = 0; ii < bonds.size(); ii++) {
            bonds[ii]->_SSSRVisited = false;
            bonds[ii]->_SSSRParent = NULL;
            bonds[ii]->_SSSRParentAtStart = true;
        }

        sketcherMinimizerBond* bond = bonds[i];

        std::queue<sketcherMinimizerBond*> q;
        bond->_SSSRVisited = true;
        q.push(bond);
        bool closedRing = false;
        while (!q.empty() && !closedRing) {
            sketcherMinimizerBond* lastBond = q.front();
            q.pop();

            sketcherMinimizerAtom* pivotAtom = lastBond->endAtom;
            if (!lastBond->_SSSRParentAtStart)
                pivotAtom = lastBond->startAtom;
            for (unsigned int j = 0; j < pivotAtom->bonds.size(); j++) {
                sketcherMinimizerBond* nextBond = pivotAtom->bonds[j];
                // sketcherMinimizerAtom *nextAtom = pivotAtom->neighbors[j];
                if (nextBond == lastBond)
                    continue;
                if (nextBond->_SSSRVisited) {
                    if (nextBond == bond) {
                        addRing(closeRing(lastBond), rings);
                        closedRing = true;
                    }

                } else {
                    if (nextBond->endAtom == pivotAtom)
                        nextBond->_SSSRParentAtStart = false;
                    nextBond->_SSSRParent = lastBond;
                    nextBond->_SSSRVisited = true;
                    q.push(nextBond);
                }
            }
        }
    }
    for (unsigned int i = 0; i < rings.size(); i++) {
        sketcherMinimizerRing* ring = rings[i];
        for (unsigned int j = 0; j < ring->_bonds.size(); j++) {
            sketcherMinimizerBond* bond = ring->_bonds[j];
            bond->rings.push_back(ring);
        }
    }
}

sketcherMinimizerRing*
sketcherMinimizerMolecule::closeRing(sketcherMinimizerBond* bond)
{
    sketcherMinimizerRing* ring = new sketcherMinimizerRing();
    sketcherMinimizerBond* lastBond = bond;
    while (lastBond) {
        ring->_bonds.push_back(lastBond);
        lastBond = lastBond->_SSSRParent;
    }
    return ring;
}

void sketcherMinimizerMolecule::addRing(
    sketcherMinimizerRing* ring, std::vector<sketcherMinimizerRing*>& rings)
{
    bool found = false;
    for (unsigned int i = 0; i < rings.size(); i++) {
        if (rings[i]->sameAs(ring)) {
            found = true;
            break;
        }
    }
    if (!found) {
        rings.push_back(ring);
    } else {
        delete ring;
    }
}
