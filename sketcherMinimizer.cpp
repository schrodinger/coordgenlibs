/*
 *  sketcherMinimizer.cpp
 *  2d_sketcher
 *
 *  Created by Nicola Zonta on 13/04/2010.
 *  Copyright Schrodinger, LLC. All rights reserved
 *
 */

#include "sketcherMinimizer.h"
#include "sketcherMinimizerMaths.h"

#include "sketcherMinimizerStretchInteraction.h"
#include "sketcherMinimizerBendInteraction.h"
#include "sketcherMinimizerClashInteraction.h"
#include <queue>
#include <stack>
#include <algorithm>
#include "CoordgenFragmenter.h"
#include "CoordgenMacrocycleBuilder.h"
#include <maeparser/Reader.hpp>
#include <fstream>

using namespace std;


#define RESIDUE_CLASH_DISTANCE_SQUARED 2.0 * 2.0

#ifndef M_PI
#define M_PI 3.1415926535897931
#endif

const int bondLength = BONDLENGTH;

static const unsigned int MINIMUM_LIGAND_ATOMS = 8;
static const float SCORE_MULTIPLIER_FOR_DOUBLE_BONDS = 0.82;
static const float SCORE_MULTIPLIER_FOR_SINGLE_BONDED_HETEROATOMS = 0.9;
static const float SCORE_MULTIPLIER_FOR_FRAGMENTS = 0.1;
const int MAX_NUMBER_OF_RINGS = 40;

sketcherMinimizer::sketcherMinimizer(float precision)
{
    m_fragmentBuilder.m_evenAngles = false;
    m_minimizer.m_evenAngles = false;
    m_minimizer.setPrecision(precision);
    m_fragmentBuilder.setPrecision(precision);
}

sketcherMinimizer::~sketcherMinimizer()
{
    clear();
}

void sketcherMinimizer::setScoreResidueInteractions(bool b)
{
    m_minimizer.m_scoreResidueInteractions = b;
}

void sketcherMinimizer::canonicalOrdering(sketcherMinimizerMolecule* minMol)
{
    vector<int> scores;
    for (unsigned int i = 0; i < minMol->_atoms.size(); i++) {
        minMol->_atoms[i]->_generalUseN = i;
    }

    sketcherMinimizer::morganScores(minMol->_atoms, minMol->_bonds, scores);

    if (scores.size() == minMol->_atoms.size()) {

        for (unsigned int i = 0; i < scores.size(); i++) {
            scores[i] *= 100;
            scores[i] += minMol->_atoms[i]->atomicNumber;
        }
        foreach (sketcherMinimizerAtom* at, minMol->_atoms) {
            at->neighbors.clear();
            at->bonds.clear();
        }
        foreach (sketcherMinimizerBond* bo, minMol->_bonds) {
            bo->startAtom->neighbors.push_back(bo->endAtom);
            bo->endAtom->neighbors.push_back(bo->startAtom);
            bo->startAtom->bonds.push_back(bo);
            bo->endAtom->bonds.push_back(bo);
        }

        vector<sketcherMinimizerAtom*> newAtoms;
        vector<sketcherMinimizerBond*> newBonds;

        for (unsigned int i = 0; i < minMol->_atoms.size(); i++) {
            minMol->_atoms[i]->_generalUseN = i;
            minMol->_atoms[i]->_generalUseVisited = false;
        }
        for (unsigned int i = 0; i < minMol->_bonds.size(); i++) {
            minMol->_bonds[i]->_SSSRVisited = false;
        }
        bool found = true;
        do {
            int scoreMaxI = -1;
            for (unsigned int i = 0; i < scores.size(); i++) {
                if (minMol->_atoms[i]->_generalUseVisited)
                    continue;
                if (scoreMaxI == -1)
                    scoreMaxI = i;
                else if (scores[i] > scores[scoreMaxI])
                    scoreMaxI = i;
            }
            if (scoreMaxI > -1) {
                queue<sketcherMinimizerAtom*> q;
                q.push(minMol->_atoms[scoreMaxI]);
                minMol->_atoms[scoreMaxI]->_generalUseVisited = true;

                while (q.size()) {
                    sketcherMinimizerAtom* at = q.front();
                    newAtoms.push_back(at);
                    q.pop();
                    int neighI = -1;
                    do {
                        neighI = -1;
                        for (unsigned int i = 0; i < at->neighbors.size();
                             i++) {
                            if (at->bonds[i]->_SSSRVisited)
                                continue;
                            else {
                                if (neighI == -1)
                                    neighI = i;
                                else {
                                    if (scores[at->neighbors[neighI]
                                                   ->_generalUseN] <
                                        scores[at->neighbors[i]
                                                   ->_generalUseN]) {
                                        neighI = i;
                                    }
                                }
                            }
                        }
                        if (neighI > -1) {
                            if (!at->neighbors[neighI]->_generalUseVisited) {
                                at->neighbors[neighI]->_generalUseVisited =
                                    true;
                                q.push(at->neighbors[neighI]);
                            }
                            at->bonds[neighI]->_SSSRVisited = true;
                            newBonds.push_back(at->bonds[neighI]);
                        }
                    } while (neighI > -1);
                }
            } else {
                found = false;
            }
        } while (found);

        minMol->_atoms = newAtoms;
        minMol->_bonds = newBonds;
    }
}

void sketcherMinimizer::initialize(
    sketcherMinimizerMolecule* minMol) // min mol is split into molecules if
                                       // needed and then added to the minimizer
{
    clear();
    _referenceAtoms = minMol->_atoms;
    _referenceBonds = minMol->_bonds;

    for (unsigned int bb = 0; bb < minMol->_bonds.size(); bb++) {
        if (minMol->_bonds[bb]->skip)
            continue;
        if (minMol->_bonds[bb]->bondOrder == 1 ||
            minMol->_bonds[bb]->bondOrder == 2) {
            if (sketcherMinimizerAtom::isMetal(
                    minMol->_bonds[bb]->startAtom->atomicNumber) ||
                sketcherMinimizerAtom::isMetal(
                    minMol->_bonds[bb]->endAtom->atomicNumber)) {
                minMol->_bonds[bb]->bondOrder = 0;
            }
        }
    }

    for (unsigned int bb = 0; bb < minMol->_bonds.size(); bb++) {
        if (minMol->_bonds[bb]->skip)
            continue;
        if (minMol->_bonds[bb]->bondOrder == 0) {
            m_proximityRelations.push_back(minMol->_bonds[bb]);
        } else if (minMol->_bonds[bb]->isResidueInteraction()) {
            if (!minMol->_bonds[bb]->startAtom->isResidue() &&
                !minMol->_bonds[bb]->endAtom->isResidue())
                m_proximityRelations.push_back(minMol->_bonds[bb]);
        }
    }
    for (unsigned int bb = 0; bb < m_extraBonds.size(); bb++) {
        if (m_extraBonds[bb]->skip)
            continue;
        if (m_extraBonds[bb]->bondOrder == 0) {
            m_proximityRelations.push_back(m_extraBonds[bb]);
        } else if (m_extraBonds[bb]->isResidueInteraction()) {
            if (!m_extraBonds[bb]->startAtom->isResidue() &&
                !m_extraBonds[bb]->endAtom->isResidue())
                m_proximityRelations.push_back(m_extraBonds[bb]);
        }
    }

    for (unsigned int bb = 0; bb < minMol->_bonds.size(); bb++) {
        if (minMol->_bonds[bb]->skip || minMol->_bonds[bb]->bondOrder == 0) {
            minMol->_bonds.erase(minMol->_bonds.begin() + bb);

            bb--;
        } else if (minMol->_bonds[bb]->startAtom->hidden ||
                   minMol->_bonds[bb]->endAtom->hidden) {
            minMol->_bonds.erase(minMol->_bonds.begin() + bb);

            bb--;
        }
    }

    for (unsigned int aa = 0; aa < minMol->_atoms.size(); aa++) {
        if (minMol->_atoms[aa]->hidden) {
            minMol->_atoms.erase(minMol->_atoms.begin() + aa);
            aa--;
        }
    }

    canonicalOrdering(minMol); // order atoms and bonds using morgan indices to
                               // make the result input independent

    foreach (sketcherMinimizerAtom* a, minMol->_atoms) {
        if (!a->hidden) {
            _atoms.push_back(a);
        }
        if (a->isResidue())
            _residues.push_back(static_cast<sketcherMinimizerResidue*>(a));
    }

    foreach (sketcherMinimizerBond* b, minMol->_bonds) {
        if (!b->startAtom->hidden && !b->endAtom->hidden) {
            _bonds.push_back(b);
        }
        if (b->isResidueInteraction())
            _residueInteractions.push_back(
                static_cast<sketcherMinimizerResidueInteraction*>(b));
    }

    minMol->forceUpdateStruct(minMol->_atoms, minMol->_bonds, minMol->_rings);
    splitIntoMolecules(minMol, _molecules);

    foreach (sketcherMinimizerBond* b, m_proximityRelations) {
        b->startAtom->molecule->m_proximityRelations.push_back(b);
        if (b->endAtom != b->startAtom)
            b->endAtom->molecule->m_proximityRelations.push_back(b);
    }

    flagCrossAtoms();

    m_minimizer._atoms = _atoms;
    m_minimizer._bonds = _bonds;
    m_minimizer._molecules = _molecules;
    m_minimizer._residues = _residues;
    m_minimizer._residueInteractions = _residueInteractions;
}

bool sketcherMinimizer::structurePassSanityCheck() const
{
    if (!_atoms.size())
        return false;
    for (auto molecule : _molecules) {
        if (molecule->_rings.size() > MAX_NUMBER_OF_RINGS) {
            return false;
        }
    }
    return true;
}

bool sketcherMinimizer::runGenerateCoordinates()
{
    bool cleanPose = true;
    if (structurePassSanityCheck()) {
        findFragments();
        m_minimizer.buildFromFragments(true);
        cleanPose = m_minimizer.avoidClashes();
        bestRotation();
        maybeFlip();
        arrangeMultipleMolecules();
        writeStereoChemistry();
    }
    return cleanPose;
}


void sketcherMinimizer::flagCrossAtoms()
{
    foreach (sketcherMinimizerAtom* at, _atoms)
        if (at->atomicNumber == 16 || at->atomicNumber == 15)
            at->crossLayout = true;

    foreach (sketcherMinimizerAtom* at, _atoms) {
        if (at->crossLayout)
            continue;
        int cross = 0;
        foreach (sketcherMinimizerAtom* n, at->neighbors) {
            if (n->neighbors.size() > 3) {
                cross++;
            }
        }
        if (cross > 2) {
            at->crossLayout = true;
        }
    }
}

/*
void sketcherMinimizer::exportCoordinates(ChmMol& molecule)
{

    foreach (sketcherMinimizerAtom* a, _referenceAtoms) {
        a->setCoordinates(a->coordinates); // round the coordinates

        // set NaN coordinates to (0);
        if (a->coordinates.x() != a->coordinates.x()) {
            cerr << "coordgen warning: NaN coordinates" << endl;
            a->coordinates.setX(0.f);
        }
        if (a->coordinates.y() != a->coordinates.y()) {
            cerr << "coordgen warning: NaN coordinates" << endl;
            a->coordinates.setY(0.f);
        }
    }

    ChmAtoms as = molecule.getAtoms(true);
    unsigned int i = 0;

    while (as.hasNext()) {

        ChmAtom& chma = as.next();
        //   if (chma.isDummy ()) continue;
        if (i < _referenceAtoms.size() && !_referenceAtoms[i]->hidden) {
            sketcherMinimizerPointF& p = _referenceAtoms[i]->coordinates;
            ChmPoint pt(p.x() / 35.f, p.y() / 35.f, 0.0f);
            if (pt.isZero()) {
                // apply gaussian random deviate to coordinates placed at origin
                ChmPoint offset((float) m_random.nextGaussian(0.0, 1e-20),
                                (float) m_random.nextGaussian(0.0, 1e-20), 0);
                pt += offset;
            }
            chma.setCoords(pt);
        }
        i++;
    }

    unsigned int j = 0;
    ChmBonds bs = molecule.getBonds(true);
    while (bs.hasNext()) {
        ChmBond& chmb = bs.next();
        if (chmb.atom1().isDummy() && !chmb.atom1().isWildcard())
            continue;
        if (chmb.atom2().isDummy() && !chmb.atom2().isWildcard())
            continue;

        if (chmb.isHidden()) {
            chmb.setBondDirectionHint(ChmBond::directionNone);
            continue; // do not increment j
        }
        const bool stereoBonds = _referenceBonds[j]->hasStereochemistryDisplay;
        if (stereoBonds) {
            if (_referenceBonds[j]->isWedge) {
                chmb.setBondDirectionHint(ChmBond::directionUp);
            } else {
                chmb.setBondDirectionHint(ChmBond::directionDown);
            }
        } else {
            chmb.setBondDirectionHint(ChmBond::directionNone);
        }
        j++;
    }
}
*/

void sketcherMinimizer::clear()
{
    for (unsigned int i = 0; i < _referenceAtoms.size(); i++) {
        delete _referenceAtoms[i];
    }
    _referenceAtoms.clear();
    _residues.clear();

    for (unsigned int i = 0; i < _referenceBonds.size(); i++) {
        delete _referenceBonds[i];
    }

    _referenceBonds.clear();

    for (unsigned int i = 0; i < m_extraBonds.size(); i++) {
        delete m_extraBonds[i];
    }

    m_extraBonds.clear();

    for (unsigned int i = 0; i < _fragments.size(); i++) {
        delete _fragments[i];
    }

    _fragments.clear();

    for (unsigned int i = 0; i < _molecules.size(); i++) {

        delete _molecules[i];
    }

    _molecules.clear();
}


/*
void sketcherMinimizer::initializeFromMolecule(ChmMol& mol)
{
    sketcherMinimizerMolecule* minMol = new sketcherMinimizerMolecule;
    ChmAtoms as = mol.getAtoms(true);
    unsigned int i = 0;
    while (as.hasNext()) {
        ChmAtom& chma = as.next();
        sketcherMinimizerAtom* min_at = new sketcherMinimizerAtom;
        min_at->_generalUseN = i;
        min_at->_generalUseN2 = i; // generalUseN will be changed by initialise,
                                   // we need a second counter to assign atom
                                   // chirality
        min_at->m_chmN = chma.getMolIndex();

        i++;

        min_at->charge = chma.getFormalCharge();
        min_at->atomicNumber = chma.getAtomicNumber();
        min_at->fixed = chma.isFixed();

        min_at->constrained = chma.isConstrained();
        min_at->coordinates =
            sketcherMinimizerPointF(chma.getX() * 35.f, chma.getY() * 35.f);
        min_at->templateCoordinates = min_at->coordinates;

        min_at->hidden = chma.isHidden();
        if (chma.isDummy() && !chma.isWildcard())
            min_at->hidden = true;

        min_at->molecule = minMol;
        minMol->_atoms.push_back(min_at);
    }

    ChmBonds bs = mol.getBonds(true);

    while (bs.hasNext()) {
        ChmBond& chmb = bs.next();
        // do not create sketcher bond to any hidden atoms
        // this means that bond vectors will be out of sequence
        if (chmb.atom1().isHidden() || chmb.atom2().isHidden())
            continue;
        if ((chmb.atom1().isDummy() && !chmb.atom1().isWildcard()) ||
            (chmb.atom2().isDummy() && !chmb.atom2().isWildcard()))
            continue;
        sketcherMinimizerBond* min_bo = new sketcherMinimizerBond;
        min_bo->startAtom = minMol->_atoms[chmb.atom1().getMolIndex()];
        min_bo->endAtom = minMol->_atoms[chmb.atom2().getMolIndex()];
        int order = chmb.getOrder();
        if (chmb.isAromatic())
            order -= 4;
        min_bo->bondOrder = order;
        minMol->_bonds.push_back(min_bo);
    }

    vector<pair<ATOM_INDEX_TYPE, ATOM_INDEX_TYPE>> zobs =
        mol.getZeroOrderBonds();
    const unsigned int zbsz = (int) zobs.size();
    unsigned int atomsize = (int) as.size();
    if (minMol->_atoms.size() < atomsize)
        atomsize = (int) minMol->_atoms.size();
    for (unsigned int i = 0; i < zbsz; i++) {
        pair<ATOM_INDEX_TYPE, ATOM_INDEX_TYPE>& pr = zobs[i];
        const unsigned int i1 = pr.first;
        const unsigned int i2 = pr.second;
        if (i1 >= atomsize || i2 >= atomsize)
            continue;
        ChmAtom& atom1 = as[i1];
        ChmAtom& atom2 = as[i2];
        if (atom1.isHidden() || atom2.isHidden())
            continue;
        if ((atom1.isDummy() && !atom1.isWildcard()) ||
            (atom2.isDummy() && !atom2.isWildcard()))
            continue;
        sketcherMinimizerBond* min_bo = new sketcherMinimizerBond;
        min_bo->startAtom = minMol->_atoms[i1];
        min_bo->endAtom = minMol->_atoms[i2];
        min_bo->bondOrder = 0;
        minMol->_bonds.push_back(min_bo);
    }

    sketcherMinimizerMolecule::assignBondsAndNeighbors(minMol->_atoms,
                                                       minMol->_bonds);

    foreach (sketcherMinimizerBond* b, minMol->_bonds) {
        if (b->bondOrder == 2) {

            sketcherMinimizerAtom* startA = b->startAtomCIPFirstNeighbor();
            sketcherMinimizerAtom* endA = b->endAtomCIPFirstNeighbor();
            if (startA && endA) {
                ChmBond chb =
                    mol.getBond(mol.getAtom(b->startAtom->_generalUseN2),
                                mol.getAtom(b->endAtom->_generalUseN2));

                if (chb.getStereo() == BondCis) {
                    int i1 = -1, i2 = -1;
                    chb.getCis(i1, i2);
                    bool isZ = true;
                    if (!(i1 == startA->_generalUseN2 ||
                          i1 == endA->_generalUseN2))
                        isZ = !isZ;
                    if (!(i2 == startA->_generalUseN2 ||
                          i2 == endA->_generalUseN2))
                        isZ = !isZ;
                    b->isZ = isZ;
                    //   cerr << "translating it to "<<isZ<<endl;
                } else if (chb.getStereo() == BondTrans) {
                    int i1 = -1, i2 = -1;
                    chb.getTrans(i1, i2);
                    bool isZ = false;
                    if (!(i1 == startA->_generalUseN2 ||
                          i1 == endA->_generalUseN2))
                        isZ = !isZ;
                    if (!(i2 == startA->_generalUseN2 ||
                          i2 == endA->_generalUseN2))
                        isZ = !isZ;
                    b->isZ = isZ;
                } else {
                    b->m_ignoreZE = true;
                }
            }
        }
    }
    initialize(minMol);
    for (unsigned int j = 0; j < minMol->_atoms.size(); j++) {
        sketcherMinimizerAtom* a = minMol->_atoms[j];

        ChmAtom atom = mol.getAtom(a->_generalUseN2);
        ChmChiralityInfo info = atom.getChiralityInfo();
        a->setStereochemistryFromChmChiralityInfo(info);
    }
}
*/

void sketcherMinimizer::splitIntoMolecules(
    sketcherMinimizerMolecule* mol, vector<sketcherMinimizerMolecule*>& mols)
{

    if (!mol->_atoms.size()) {
        mols.push_back(mol);
        return;
    }
    foreach (sketcherMinimizerAtom* a, mol->_atoms)
        a->_generalUseVisited = false;
    queue<sketcherMinimizerAtom*> q;
    q.push(mol->_atoms[0]);
    foreach (sketcherMinimizerAtom* a, mol->_atoms) {
        if (!a->hidden) {
            q.push(a);
            break;
        }
    }
    while (q.size()) {
        sketcherMinimizerAtom* a = q.front();
        q.pop();
        a->_generalUseVisited = true;
        foreach (sketcherMinimizerAtom* n, a->neighbors) {
            if (!n->_generalUseVisited && !n->hidden)
                q.push(n);
        }
    }
    vector<sketcherMinimizerAtom*> newAtoms;

    foreach (sketcherMinimizerAtom* a, mol->_atoms) {
        if (!a->_generalUseVisited && !a->hidden) {
            newAtoms.push_back(a);
        }
    }
    if (!newAtoms.size()) {
        mols.push_back(mol);
        foreach (sketcherMinimizerMolecule* m, mols) {
            foreach (sketcherMinimizerAtom* a, m->_atoms) {
                a->_generalUseVisited = false;
            }
        }

    } else {
        sketcherMinimizerMolecule* newMol = new sketcherMinimizerMolecule;
        for (unsigned int i = 0; i < mol->_rings.size(); i++) {

            if (!mol->_rings[i]->_atoms[0]->_generalUseVisited) {
                newMol->_rings.push_back(mol->_rings[i]);
                mol->_rings.erase(mol->_rings.begin() + i);
                i--;
            }
        }
        for (unsigned int i = 0; i < mol->_bonds.size(); i++) {
            if (!mol->_bonds[i]->startAtom->_generalUseVisited) {
                newMol->_bonds.push_back(mol->_bonds[i]);
                mol->_bonds.erase(mol->_bonds.begin() + i);
                i--;
            }
        }

        for (unsigned int i = 0; i < mol->_atoms.size(); i++) {
            if (!mol->_atoms[i]->_generalUseVisited) {
                mol->_atoms[i]->molecule = newMol;
                newMol->_atoms.push_back(mol->_atoms[i]);
                mol->_atoms.erase(mol->_atoms.begin() + i);
                i--;
            }
        }
        mols.push_back(mol);
        splitIntoMolecules(newMol, mols);
    }
}

sketcherMinimizerRing*
sketcherMinimizer::sameRing(const sketcherMinimizerAtom* at1,
                            const sketcherMinimizerAtom* at2)
{
    return sketcherMinimizerAtom::shareARing(at1, at2);
}

sketcherMinimizerRing*
sketcherMinimizer::sameRing(const sketcherMinimizerAtom* at1,
                            const sketcherMinimizerAtom* at2,
                            const sketcherMinimizerAtom* at3)
{
    if (!at1->rings.size())
        return NULL;
    if (!at2->rings.size())
        return NULL;
    if (!at3->rings.size())
        return NULL;
    sketcherMinimizerRing* r = 0;
    foreach (sketcherMinimizerRing* ring, at1->rings) {
        if (ring->isMacrocycle())
            continue;
        foreach (sketcherMinimizerRing* ring2, at2->rings) {
            if (ring != ring2)
                continue;
            foreach (sketcherMinimizerRing* ring3, at3->rings) {
                if (ring3 == ring2) {
                    if (!r)
                        r = ring2;
                    else if (ring2->_atoms.size() < r->_atoms.size())
                        r = ring2;
                }
            }
        }
    }
    foreach (sketcherMinimizerRing* ring, at1->rings) {
        foreach (sketcherMinimizerRing* ring2, at2->rings) {
            if (ring != ring2)
                continue;
            foreach (sketcherMinimizerRing* ring3, at3->rings) {
                if (ring3 == ring2) {
                    if (!r)
                        r = ring2;
                    else if (ring2->_atoms.size() < r->_atoms.size())
                        r = ring2;
                }
            }
        }
    }
    return r;
}

void sketcherMinimizer::writeStereoChemistry()
{
    foreach (sketcherMinimizerAtom* a, _atoms) {
        if (a->hasStereochemistrySet)
            a->writeStereoChemistry();
    }
    assignPseudoZ();
}

void sketcherMinimizer::assignPseudoZ()
{
    foreach (sketcherMinimizerMolecule* mol, _molecules) {
        foreach (sketcherMinimizerAtom* a, mol->_atoms) {
            a->_generalUseVisited = false;
        }
        sketcherMinimizerAtom* lastAtom = NULL;
        bool finished = false;
        while (!finished) {
            lastAtom = NULL;
            foreach (sketcherMinimizerAtom* a, mol->_atoms) {
                if (!a->_generalUseVisited) {
                    lastAtom = a;
                    break;
                }
            }
            if (lastAtom) {
                queue<sketcherMinimizerAtom*> q;
                q.push(lastAtom);
                while (q.size()) {
                    lastAtom = q.front();
                    q.pop();
                    lastAtom->_generalUseVisited = true;
                    for (unsigned int i = 0; i < lastAtom->neighbors.size();
                         i++) {
                        if (lastAtom->neighbors[i]->_generalUseVisited)
                            continue;
                        float Z = lastAtom->m_pseudoZ;
                        sketcherMinimizerBond* b = lastAtom->bonds[i];
                        if (b->hasStereochemistryDisplay) {
                            if (b->isWedge) {
                                if ((b->startAtom == lastAtom &&
                                     b->isReversed == false) ||
                                    (b->endAtom == lastAtom &&
                                     b->isReversed == true)) {
                                    Z += 1.f;
                                } else if ((b->startAtom == lastAtom &&
                                            b->isReversed == true) ||
                                           (b->endAtom == lastAtom &&
                                            b->isReversed == false)) {
                                    Z -= 1.f;
                                }

                            } else {
                                if ((b->startAtom == lastAtom &&
                                     b->isReversed == false) ||
                                    (b->endAtom == lastAtom &&
                                     b->isReversed == true)) {
                                    Z -= 1.f;
                                } else if ((b->startAtom == lastAtom &&
                                            b->isReversed == true) ||
                                           (b->endAtom == lastAtom &&
                                            b->isReversed == false)) {
                                    Z += 1.f;
                                }
                            }
                        }
                        lastAtom->neighbors[i]->m_pseudoZ = Z;
                        q.push(lastAtom->neighbors[i]);
                    }
                }
            } else
                finished = true;
        }
    }
}

void sketcherMinimizer::maybeFlipPeptides(
    std::vector<sketcherMinimizerAtom*> atoms, float& scoreX)
{
    auto chetoCs = m_minimizer.getChetoCs(atoms);
    auto aminoNs = m_minimizer.getAminoNs(atoms);
    auto alphaCs = m_minimizer.getAlphaCs(atoms, chetoCs, aminoNs);
    for (auto alphaC : alphaCs) {
        sketcherMinimizerAtom* aminoN = nullptr;
        sketcherMinimizerAtom* chetoC = nullptr;
        for (auto neighbor : alphaC->neighbors) {
            if (aminoNs.find(neighbor) != aminoNs.end()) {
                aminoN = neighbor;
            } else if (chetoCs.find(neighbor) != chetoCs.end()) {
                chetoC = neighbor;
            }
        }
        if (aminoN && chetoC) {
            auto direction = aminoN->coordinates - chetoC->coordinates;
            const float PEPTIDE_SCORE = 100.f;
            if (direction.x() > 0) {
                scoreX -= PEPTIDE_SCORE;
            } else {
                scoreX += PEPTIDE_SCORE;
            }
        }
    }
}

void sketcherMinimizer::maybeFlip()
{
    foreach (sketcherMinimizerMolecule* mol, _molecules) {
        if (mol->hasFixedFragments)
            continue;
        if (mol->hasConstrainedFragments)
            continue;
        if (mol->_atoms.size() < 2)
            continue;

        float scoreY = 0.f, scoreX = 0.f;
        maybeFlipPeptides(mol->getAtoms(), scoreX);
        sketcherMinimizerPointF cent(0.f, 0.f);
        foreach (sketcherMinimizerAtom* a, mol->_atoms)
            cent += a->coordinates;
        if (mol->_atoms.size())
            cent /= mol->_atoms.size();

        foreach (sketcherMinimizerFragment* f, mol->_fragments) {
            vector<sketcherMinimizerRing*> rings = f->getRings();
            if (rings.size() == 2) {
                sketcherMinimizerRing* r1 = rings[0];
                sketcherMinimizerRing* r2 = rings[1];
                sketcherMinimizerPointF c1 = r1->findCenter();
                sketcherMinimizerPointF c2 = r2->findCenter();
                if (c1.x() - c2.x() > SKETCHER_EPSILON) {
                    sketcherMinimizerPointF swapP = c2;
                    c2 = c1;
                    c1 = swapP;
                    sketcherMinimizerRing* swapR;
                    swapR = r2;
                    r2 = r1;
                    r1 = swapR;
                }
                if (r2->isBenzene() && !r1->isBenzene()) {
                    scoreX -= 20;
                } else if (r1->isBenzene() && !r2->isBenzene()) {
                    scoreX += 20;
                }
            }

            if (rings.size() > 3) {

                sketcherMinimizerPointF center(0.f, 0.f);
                sketcherMinimizerPointF weightedCenter(0.f, 0.f);

                int totalN = 0;
                int totalRings = 0;
                foreach (sketcherMinimizerRing* r, rings) {
                    if (r->_atoms.size() < 4)
                        continue;

                    sketcherMinimizerPointF c = r->findCenter();
                    center += c;
                    weightedCenter += c * r->_atoms.size();
                    totalN += r->_atoms.size();
                    totalRings++;
                }
                if (totalRings && totalN) {
                    center /= totalRings;
                    weightedCenter /= totalN;
                }
                if (weightedCenter.y() - center.y() < -SKETCHER_EPSILON) {
                    scoreY += 50.f;
                } else if (weightedCenter.y() - center.y() > SKETCHER_EPSILON) {
                    scoreY -= 50.f;
                }

                if (weightedCenter.x() - center.x() < -SKETCHER_EPSILON) {
                    scoreX += 50.f;
                } else if (weightedCenter.x() - center.x() > SKETCHER_EPSILON) {
                    scoreX -= 50.f;
                }
            }
        }

        float minx = 9999.f, miny = 9999.f, maxx = -9999.f, maxy = -9999.f;

        foreach (sketcherMinimizerAtom* a, mol->_atoms) {
            float x = a->coordinates.x();
            float y = a->coordinates.y();
            if (x < minx)
                minx = x;
            else if (x > maxx)
                maxx = x;
            if (y < miny)
                miny = y;
            else if (y > maxy)
                maxy = y;
        }

        float meanx = (maxx + minx) * 0.5;
        float meany = (maxy + miny) * 0.5;

        if (meanx - cent.x() > SKETCHER_EPSILON)
            scoreX -= 0.5;
        else if (meanx - cent.x() < -SKETCHER_EPSILON)
            scoreX += 0.5;
        if (meany - cent.y() > SKETCHER_EPSILON)
            scoreY += 0.5;
        else if (meany - cent.y() < -SKETCHER_EPSILON)
            scoreY -= 0.5;

        foreach (sketcherMinimizerBond* b, mol->_bonds) {
            if (b->bondOrder == 2) {

                if (b->startAtom->neighbors.size() == 1 &&
                    b->endAtom->neighbors.size() > 1) {

                    float diff = b->startAtom->coordinates.y() -
                                 b->endAtom->coordinates.y();
                    if (diff > SKETCHER_EPSILON)
                        scoreY += 1;
                    else if (diff < -SKETCHER_EPSILON)
                        scoreY -= 1;
                } else if (b->endAtom->neighbors.size() == 1 &&
                           b->startAtom->neighbors.size() > 1) {

                    float diff = b->endAtom->coordinates.y() -
                                 b->startAtom->coordinates.y();
                    if (diff > SKETCHER_EPSILON)
                        scoreY += 1;
                    else if (diff < -SKETCHER_EPSILON)
                        scoreY -= 1;
                }
            }
        }

        if (0.f > scoreY) {
            foreach (sketcherMinimizerAtom* a, mol->_atoms) {
                a->coordinates.setY(-a->coordinates.y());
            }
        }
        if (0.f > scoreX) {

            foreach (sketcherMinimizerAtom* a, mol->_atoms) {
                a->coordinates.setX(-a->coordinates.x());
            }
        }
    }
}

void sketcherMinimizer::addToVector(float weight, float angle,
                                    vector<pair<float, float>>& angles)
{
    angle = roundToTwoDecimalDigits(angle);
    while (angle <= 0)
        angle += M_PI;
    for (unsigned int i = 0; i < angles.size(); i++) {
        if (angles[i].second < angle - SKETCHER_EPSILON) {
            if (i == angles.size() - 1) {
                angles.push_back(pair<float, float>(weight, angle));
                break;
            }
        } else if (angles[i].second - angle < SKETCHER_EPSILON &&
                   angles[i].second - angle > -SKETCHER_EPSILON) {
            angles[i].first += weight;
            break;
        } else {
            angles.insert(angles.begin() + i,
                          pair<float, float>(weight, angle));
            break;
        }
    }
    if (!angles.size())
        angles.push_back(pair<float, float>(weight, angle));
}

/*if a peptide chain is present rotate the molecule so it's horizontal*/
void sketcherMinimizer::addBestRotationInfoForPeptides(
    vector<pair<float, float>>& angles,
    std::vector<sketcherMinimizerAtom*> atoms)
{
    auto chetoCs = m_minimizer.getChetoCs(atoms);
    auto aminoNs = m_minimizer.getAminoNs(atoms);
    auto alphaCs = m_minimizer.getAlphaCs(atoms, chetoCs, aminoNs);
    for (auto alphaC : alphaCs) {
        sketcherMinimizerAtom* aminoN = nullptr;
        sketcherMinimizerAtom* chetoC = nullptr;
        for (auto neighbor : alphaC->neighbors) {
            if (aminoNs.find(neighbor) != aminoNs.end()) {
                aminoN = neighbor;
            } else if (chetoCs.find(neighbor) != chetoCs.end()) {
                chetoC = neighbor;
            }
        }
        if (aminoN && chetoC) {
            auto direction = aminoN->coordinates - chetoC->coordinates;
            float weight = 1000.f;
            float angle = atan2(-direction.y(), direction.x());
            addToVector(weight, angle, angles);
        }
    }
}

void sketcherMinimizer::bestRotation()
{
    foreach (sketcherMinimizerMolecule* mol, _molecules) {
        vector<pair<float, float>> angles;
        if (mol->hasFixedFragments)
            continue;
        if (mol->getMainFragment()) {
            if (mol->getMainFragment()->constrained)
                continue;
            if (mol->getMainFragment()->isTemplated)
                continue;
        }
        addBestRotationInfoForPeptides(angles, mol->getAtoms());
        float angle = 0.f;
        float lastAngle;
        unsigned int i = 0, j = 0;
        float weight = 1.f;
        float increment = M_PI / 6;
        foreach (sketcherMinimizerAtom* a, mol->_atoms) {
            if (a->rings.size())
                continue;
            if (a->neighbors.size() > 1) {
                for (i = 0; i < a->neighbors.size() - 1; i++) {
                    for (j = i + 1; j < a->neighbors.size(); j++) {

                        weight = 6;

                        if (a->neighbors[i]->neighbors.size() != 1)
                            weight += 2;
                        if (a->neighbors[j]->neighbors.size() != 1)
                            weight += 2;
                        if (a->neighbors[j]->atomicNumber == 6)
                            weight += 1;
                        if (a->neighbors[j]->atomicNumber == 6)
                            weight += 1;
                        if (a->neighbors[i]->charge == 0)
                            weight += 1;
                        if (a->neighbors[j]->charge == 0)
                            weight += 1;

                        sketcherMinimizerPointF p =
                            a->neighbors[i]->coordinates -
                            a->neighbors[j]->coordinates;
                        angle = atan2(-p.y(), p.x());
                        addToVector(weight, angle, angles);
                    }
                }
            }
        }
        foreach (sketcherMinimizerBond* b, mol->_bonds) {
            sketcherMinimizerPointF p =
                b->endAtom->coordinates - b->startAtom->coordinates;
            weight = 1;
            angle = atan2(-p.y(), p.x());
            angle = roundToTwoDecimalDigits(angle);

            while (angle <= 0)
                angle += M_PI;
            lastAngle = angle;
            for (unsigned int i = 0; i < 6; i++) {
                if (i == 1 || i == 5)
                    weight = 5.f;
                else if (i == 0 || i == 3)
                    weight = 1.5f;
                else
                    weight = 1.f;
                if (b->bondOrder == 2 && i == 3 &&
                    (b->startAtom->neighbors.size() == 1 ||
                     b->endAtom->neighbors.size() == 1))
                    weight += 1.5;

                if (b->startAtom->neighbors.size() == 1 &&
                    b->endAtom->neighbors.size() == 1 && i == 0)
                    weight += 10;
                addToVector(weight, lastAngle, angles);
                lastAngle += increment;
                if (lastAngle > M_PI)
                    lastAngle -= M_PI;
            }
        }

        foreach (sketcherMinimizerFragment* f, mol->_fragments) {

            vector<sketcherMinimizerRing*> rings = f->getRings();
            vector<sketcherMinimizerRing*> inPlaneRings = rings;

            int ringsN = inPlaneRings.size();
            if (ringsN == 2) {

                sketcherMinimizerRing* r1 = inPlaneRings[0];
                sketcherMinimizerRing* r2 = inPlaneRings[1];
                sketcherMinimizerPointF p = r2->findCenter() - r1->findCenter();
                p.normalize();
                angle = atan2(-p.y(), p.x());
                weight = 25.f;

                addToVector(weight, angle, angles);

            } else if (ringsN == 3) {

                sketcherMinimizerPointF r1 = inPlaneRings[0]->findCenter();
                sketcherMinimizerPointF r2 = inPlaneRings[1]->findCenter();

                foreach (sketcherMinimizerRing* r, inPlaneRings) {
                    vector<sketcherMinimizerRing*> fusedWith;
                    vector<vector<sketcherMinimizerAtom*>> fusionAtoms;

                    for (unsigned int fw = 0; fw < r->fusedWith.size(); fw++) {
                        fusedWith.push_back(r->fusedWith[fw]);
                        fusionAtoms.push_back(r->fusionAtoms[fw]);
                    }
                    if (fusedWith.size() == 2) {
                        if (fusionAtoms[0].size() == 2 &&
                            fusionAtoms[1].size() == 2) {
                            r1 = (fusionAtoms[0][0]->coordinates +
                                  fusionAtoms[0][1]->coordinates) *
                                 0.5;
                            r2 = (fusionAtoms[1][0]->coordinates +
                                  fusionAtoms[1][1]->coordinates) *
                                 0.5;
                            break;
                        }
                    }
                }

                sketcherMinimizerPointF p = r2 - r1;
                angle = atan2(-p.y(), p.x());
                weight = 50.f;
                addToVector(weight, angle, angles);

            } else {
                vector<sketcherMinimizerRing*> rings;
                foreach (sketcherMinimizerRing* r, inPlaneRings) {
                    if (r->_atoms.size() == 6) {
                        rings.push_back(r);
                    }
                }
                foreach (sketcherMinimizerRing* r, rings) {
                    for (unsigned int i = 0; i < r->fusionAtoms.size(); i++) {
                        vector<sketcherMinimizerAtom*> fusionAts =
                            r->fusionAtoms[i];
                        if (fusionAts.size() == 2) {
                            sketcherMinimizerPointF p =
                                fusionAts[0]->coordinates -
                                fusionAts[1]->coordinates;
                            //  if (p.x () != p.x () || p.y () != p.y ()) p =
                            //  sketcherMinimizerPointF (50.f, 0.f);
                            sketcherMinimizerPointF rotatedP(p.y(), p.x());
                            angle = atan2(-p.y(), p.x()) - M_PI * 0.5;
                            weight = 25.f;
                            addToVector(weight, angle, angles);
                        }
                    }
                }
            }
        }
        if (angles.size() > 1) {
            if (angles[angles.size() - 1].second - angles[0].second >=
                M_PI - 2 * SKETCHER_EPSILON) {
                angles[0].first += angles[angles.size() - 1].first;
                angles.erase(angles.begin() + angles.size() - 1);
            }
        }

        if (angles.size()) {
            int bestI = 0;
            for (i = 0; i < angles.size(); i++) {
                if (angles[i].first > angles[bestI].first)
                    bestI = i;
            }
            float s = -sin(angles[bestI].second);
            float c = cos(angles[bestI].second);
            sketcherMinimizerPointF center(0.f, 0.f);
            foreach (sketcherMinimizerAtom* at, mol->_atoms)
                center += at->coordinates;
            if (mol->_atoms.size())
                center /= mol->_atoms.size();

            foreach (sketcherMinimizerAtom* at, mol->_atoms) {
                sketcherMinimizerPointF v = at->coordinates - center;
                v.rotate(s, c);
                at->setCoordinates(center + v);
            }
        }
    }
}

void sketcherMinimizer::findFragments()
{

    assert(_molecules.size());
    foreach (sketcherMinimizerMolecule* mol, _molecules) {
        CoordgenFragmenter::splitIntoFragments(mol);
        if (!mol->_fragments.size())
            continue;
        vector<sketcherMinimizerFragment*> fragments = mol->_fragments;
        _fragments.reserve(_fragments.size() + fragments.size());
        _fragments.insert(_fragments.end(), fragments.begin(), fragments.end());
        _independentFragments.push_back(mol->getMainFragment());
    }

    m_minimizer._fragments = _fragments;

    initializeFragments();
}

void sketcherMinimizer::placeResiduesProteinOnlyModeCircleStyle(
    std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains)
{
    int totalResiduesNumber = _residues.size() + chains.size();

    float angle = 2.f * M_PI / totalResiduesNumber;
    const float residueRadius = 30;
    const float circumference = totalResiduesNumber * residueRadius * 2;
    const float radius = circumference * 0.5 / M_PI;
    int i = 0;
    for (auto chain : chains) {
        ++i; // gap between chains
        auto residues = chain.second;
        sort(residues.begin(), residues.end(),
             [](const sketcherMinimizerResidue* firstRes,
                const sketcherMinimizerResidue* secondRes) {
                 int firstN = firstRes->resnum;
                 int secondN = secondRes->resnum;
                 return firstN < secondN;
             });

        for (auto res : residues) {
            sketcherMinimizerPointF p(radius, 0);
            // place residues in a circle
            p.rotate(sin(angle * i), cos(angle * i));
            res->coordinates = p;
            res->coordinatesSet = true;
            res->molecule->isPlaced = true;
            ++i;
        }
    }
}

std::map<std::string, sketcherMinimizerPointF>
sketcherMinimizer::computeChainsStartingPositionsMetaMol(
    std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains)
{
    map<std::string, sketcherMinimizerAtom*> molMap;

    sketcherMinimizerMolecule* metaMol = new sketcherMinimizerMolecule;
    for (auto chainPair : chains) {
        sketcherMinimizerAtom* a = new sketcherMinimizerAtom;
        a->molecule = metaMol;
        metaMol->_atoms.push_back(a);
        molMap[chainPair.first] = a;
    }

    for (auto chainPair : chains) {
        for (auto residue : chainPair.second) {
            for (auto interaction : residue->residueInteractions) {
                if (interaction->startAtom->isResidue() &&
                    interaction->endAtom->isResidue()) {
                    sketcherMinimizerResidue* r1 =
                        static_cast<sketcherMinimizerResidue*>(
                            interaction->startAtom);
                    sketcherMinimizerResidue* r2 =
                        static_cast<sketcherMinimizerResidue*>(
                            interaction->endAtom);
                    if (r1->chain != r2->chain) {
                        // add a bond to the metaMol if it doesn't exist already
                        sketcherMinimizerAtom* at1 = molMap[r1->chain];
                        sketcherMinimizerAtom* at2 = molMap[r2->chain];
                        bool found = false;
                        for (sketcherMinimizerBond* b : metaMol->_bonds) {
                            if ((b->startAtom == at1 && b->endAtom == at2) ||
                                (b->startAtom == at2 && b->endAtom == at1)) {
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            sketcherMinimizerBond* newBond =
                                new sketcherMinimizerBond;
                            newBond->startAtom = at1;
                            newBond->endAtom = at2;
                            metaMol->_bonds.push_back(newBond);
                        }
                    }
                }
            }
        }
    }

    sketcherMinimizer min;
    if (metaMol->_atoms.size()) {
        min.m_fragmentBuilder.m_evenAngles = true;
        min.m_minimizer.m_evenAngles = true;

        min.initialize(metaMol);
        min.findFragments();
        min.m_minimizer.buildFromFragments(true);
        min.m_minimizer.avoidClashes();
        min.bestRotation();
        min.maybeFlip();
        min.arrangeMultipleMolecules();
    }
    std::map<std::string, sketcherMinimizerPointF> positions;
    for (auto iter : molMap) {
        positions[iter.first] = iter.second->coordinates * 10.;
    }
    return positions;
}

void sketcherMinimizer::shortenInteractions(
    std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains)
{
    for (auto chain : chains) {
        for (auto res : chain.second) {
            for (auto interaction : res->residueInteractions) {
                sketcherMinimizerPointF midPoint =
                    0.5 * (interaction->startAtom->coordinates +
                           interaction->endAtom->coordinates);
                res->coordinates += (midPoint - res->coordinates) * 0.1;
            }
        }
    }
}


std::vector<sketcherMinimizerResidue*> sketcherMinimizer::orderResiduesOfChains(
    std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains)
{
    std::vector<sketcherMinimizerResidue*> vec;
    for (auto chain : chains) {
        for (auto res : chain.second) {
            vec.push_back(res);
        }
    }
    sort(vec.begin(), vec.end(), [](const sketcherMinimizerResidue* firstRes,
                                    const sketcherMinimizerResidue* secondRes) {
        return firstRes->residueInteractions.size() >
               secondRes->residueInteractions.size();
    });
    std::set<sketcherMinimizerResidue*> visitedResidues;
    std::queue<sketcherMinimizerResidue*> residueQueue;
    std::vector<sketcherMinimizerResidue*> finalVec;
    for (auto residue : vec) {
        if (visitedResidues.find(residue) != visitedResidues.end()) {
            continue;
        }
        residueQueue.push(residue);
        visitedResidues.insert(residue);
        while (residueQueue.size()) {
            auto topResidue = residueQueue.front();
            finalVec.push_back(topResidue);
            residueQueue.pop();
            for (auto partner : topResidue->residueInteractionPartners) {
                sketcherMinimizerResidue* partnerRes =
                    static_cast<sketcherMinimizerResidue*>(partner);
                if (visitedResidues.find(partnerRes) == visitedResidues.end()) {
                    residueQueue.push(partnerRes);
                    visitedResidues.insert(partnerRes);
                }
            }
        }
    }
    return finalVec;
}

void sketcherMinimizer::placeResiduesProteinOnlyModeLIDStyle(
    std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains)
{
    auto positions = computeChainsStartingPositionsMetaMol(chains);
    sketcherMinimizerPointF p;
    for (auto chain : chains) {
        p = positions[chain.first];
        for (auto res : chain.second) {
            res->coordinates = p;
        }
    }
    shortenInteractions(chains);
    auto residues = orderResiduesOfChains(chains);
    for (auto res : residues) {
        int partnersAlreadySet = 0;
        sketcherMinimizerResidue* firstPartner = NULL;
        for (auto partner : res->residueInteractionPartners) {
            if (partner->coordinatesSet) {
                ++partnersAlreadySet;
                auto partnerResidue =
                    static_cast<sketcherMinimizerResidue*>(partner);
                if (!firstPartner && partnerResidue->chain != res->chain) {
                    firstPartner = partnerResidue;
                }
            }
        }

        /*when searching for a position for res prefer a direction perpendicular
         to the direction of interactions
         to optimize use of space*/
        sketcherMinimizerPointF direction(0, 1);
        if (firstPartner) {
            sketcherMinimizerPointF chainCenterDirections =
                positions[firstPartner->chain] - positions[res->chain];
            direction = sketcherMinimizerPointF(-chainCenterDirections.y(),
                                                chainCenterDirections.x());
            direction.normalize();
            direction *= 4.;
        }
        res->coordinates = exploreGridAround(res->coordinates, 10, 5, 0, 0,
                                             -1.f, false, res, direction);
        res->coordinatesSet = true;
        res->molecule->isPlaced = true;
    }
}

void sketcherMinimizer::placeResiduesProteinOnlyMode()
{

    std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains;
    for (auto residue : _residues) {
        string chainOfResidue = residue->chain;
        chains[chainOfResidue].push_back(residue);
    }
    placeResiduesProteinOnlyModeLIDStyle(chains);
    m_minimizer.minimizeProteinOnlyLID(chains);
}

void sketcherMinimizer::placeResiduesInCrowns()
{
    auto SSEs = groupResiduesInSSEs(_residues);
    /*sort secondary structure elements so that the most importants are placed
      first.
      prefer longer SSEs and ones that make more interactions*/
    sort(SSEs.begin(), SSEs.end(),
         [](const vector<sketcherMinimizerResidue*>& firstSSE,
            const vector<sketcherMinimizerResidue*>& secondSSE) {
             float interactionsOfFirst = 0, interactionsOfSecond = 0;
             for (auto res : firstSSE) {
                 interactionsOfFirst += res->residueInteractions.size();
             }
             for (auto res : secondSSE) {
                 interactionsOfSecond += res->residueInteractions.size();
             }
             float interactionScaling = 3.f;
             float score1 =
                 firstSSE.size() +
                 interactionScaling * interactionsOfFirst / firstSSE.size();
             float score2 =
                 secondSSE.size() +
                 interactionScaling * interactionsOfSecond / secondSSE.size();
             return score1 > score2;
         });
    bool needOtherShape = true;
    int shapeCounter = 0;

    // place residues in a crowns around the ligand. Keep expanding to further
    // away crowns
    // until all residues are placed
    while (needOtherShape) {
        vector<sketcherMinimizerPointF> shape =
            shapeAroundLigand(shapeCounter++);
        needOtherShape = fillShape(SSEs, shape, shapeCounter);
    }
}

/*place residues in SSEs in the current shape. Return false if all residues are
 * place, true otherwise*/
bool sketcherMinimizer::fillShape(
    vector<vector<sketcherMinimizerResidue*>>& SSEs,
    const vector<sketcherMinimizerPointF>& shape, int shapeN)
{
    vector<bool> penalties(shape.size(), false);
    std::set<sketcherMinimizerResidue*> outliers;
    for (auto SSE : SSEs) {
        placeSSE(SSE, shape, shapeN, penalties, outliers);
    }
    return !outliers.empty();
}

/*assign a penalty for the stretching of bonds between residues in the same
 * SSE*/
float sketcherMinimizer::scoreSSEBondStretch(
    sketcherMinimizerPointF coordinates1, sketcherMinimizerPointF coordinates2)
{
    const float stretchPenalty = 400.f;
    auto squaredLength = (coordinates2 - coordinates1).squareLength();
    return squaredLength * stretchPenalty;
}


float sketcherMinimizer::getResidueDistance(
    float startF, float increment, sketcherMinimizerResidue* resToConsider,
    vector<sketcherMinimizerResidue*> SSE)
{
    float totalF = startF;
    sketcherMinimizerResidue* lastRes = nullptr;
    for (auto res : SSE) {
        if (lastRes) {
            float result = res->resnum - lastRes->resnum;
            /*if the gap is more than 1, make the distance a bit smaller for
             * aesthetic reasons*/
            result = 1 + (result - 1) * 0.8;
            if (result < 1.f)
                result = 1.f;
            totalF += increment * result;
        }
        if (res == resToConsider)
            break;
        lastRes = res;
    }
    return totalF;
}

/*return a score for the placing on the SSE starting at startingPosition and
 * separated by increment*/
float sketcherMinimizer::scoreSSEPosition(
    vector<sketcherMinimizerResidue*> SSE,
    const vector<sketcherMinimizerPointF>& shape, int shapeN,
    vector<bool>& penalties, float startingPosition, float increment)
{
    float score = 0.f;
    sketcherMinimizerResidue* lastResidue = nullptr;
    int lastResiduePosition = 0;
    sketcherMinimizerPointF lastResidueCoordinates;
    for (auto res : SSE) {
        int index = getShapeIndex(
            shape, getResidueDistance(startingPosition, increment, res, SSE));
        auto residueCoordinates = shape.at(index);
        int residuePosition = 0;
        if (res->coordinatesSet) {
            residuePosition = -1;
            residueCoordinates = res->coordinates;
        } else {
            if (!penalties[index])
                residuePosition = 0;
            else
                residuePosition = 1;
        }
        if (residuePosition != -1) {
            score += scoreResiduePosition(index, shape, shapeN, penalties, res);
        }

        if (lastResidue && (residuePosition != lastResiduePosition))

        {

            score +=
                scoreSSEBondStretch(residueCoordinates, lastResidueCoordinates);
        }
        lastResiduePosition = residuePosition;
        lastResidueCoordinates = residueCoordinates;
        lastResidue = res;
    }
    return score;
}

void sketcherMinimizer::placeSSE(vector<sketcherMinimizerResidue*> SSE,
                                 const vector<sketcherMinimizerPointF>& shape,
                                 int shapeN, vector<bool>& penalties,
                                 set<sketcherMinimizerResidue*>& outliers,
                                 bool placeOnlyInteracting)
{
    int residuesToPlace = 0;
    for (auto res : SSE) {
        if (!res->coordinatesSet) {
            residuesToPlace++;
        }
    }
    if (residuesToPlace == 0) {
        return;
    }
    typedef pair<float, float> Solution;
    vector<pair<float, Solution>> scoredSolutions;
    /*move around the crown scoring possible solutions, varying the starting
     position
     and the separation between consecutive residues*/
    for (float f = 0.f; f < 1.f; f += 0.004f) {
        float distance = 5.f / shape.size();
        for (float increment = -1 * distance; increment <= 1 * distance;
             increment += distance) {
            if (increment == 0)
                continue;
            float score =
                scoreSSEPosition(SSE, shape, shapeN, penalties, f, increment);
            scoredSolutions.push_back(
                pair<float, Solution>(score, Solution(f, increment)));
        }
    }
    auto bestResult =
        min_element(scoredSolutions.begin(), scoredSolutions.end());
    set<sketcherMinimizerResidue*> alreadyPlaced;
    for (auto residue : SSE) {
        if (residue->coordinatesSet)
            continue; // placed in a previous crown
        float f = getResidueDistance(bestResult->second.first,
                                     bestResult->second.second, residue, SSE);
        int index = getShapeIndex(shape, f);
        bool alreadyAResidueHere = penalties.at(index);
        sketcherMinimizerPointF position = shape.at(index);
        if (alreadyAResidueHere ||
            (placeOnlyInteracting &&
             !residue->residueInteractionPartners.size())) {
            outliers.insert(residue);
        } else {
            residue->coordinates = position;
            alreadyPlaced.insert(residue);
        }
    }
    // mark the current solution to prevent other residues from being placed on
    // top of these
    markSolution(bestResult->second, SSE, shape, penalties, outliers);
    for (auto res : alreadyPlaced)
        res->coordinatesSet = true;
    for (auto res : SSE) {
        if (res->m_isWaterMap && res->m_isClashing && res->coordinatesSet &&
            res->m_closestLigandAtom != nullptr) {
            sketcherMinimizerPointF directionToLigand =
                res->m_closestLigandAtom->coordinates - res->coordinates;
            directionToLigand.normalize();
            float displacement = BONDLENGTH * 0.3;
            res->coordinates = res->m_closestLigandAtom->coordinates -
                               directionToLigand * displacement;
        }
    }
}


void sketcherMinimizer::markSolution(
    pair<float, float> solution, vector<sketcherMinimizerResidue*> SSE,
    const vector<sketcherMinimizerPointF>& shape, vector<bool>& penalties,
    set<sketcherMinimizerResidue*>& outliers)
{
    float padding = abs(solution.second) * 0.5f;
    sketcherMinimizerResidue* lastRes = nullptr;
    float lastF = 0.f;
    for (auto res : SSE) {
        if (res->coordinatesSet || (res->m_isWaterMap && res->m_isClashing) ||
            outliers.find(res) != outliers.end()) {
            lastRes = nullptr;
            lastF = 0.f;
            continue;
        }
        float f = getResidueDistance(solution.first, solution.second, res, SSE);
        int startIndex = getShapeIndex(shape, f - padding);
        int endIndex = getShapeIndex(shape, f + padding);
        for (int index = startIndex; index != endIndex;
             index = (index + 1) % shape.size()) {
            penalties.at(index) = true;
        }
        if (lastRes) {
            if (solution.second < 0) {
                std::swap(lastF, f);
            }
            int startIndex = getShapeIndex(shape, lastF);
            int endIndex = getShapeIndex(shape, f);
            for (int index = startIndex; index != endIndex;
                 index = (index + 1) % shape.size()) {
                penalties.at(index) = true;
            }
        }
        lastRes = res;
        lastF = f;
    }
}

int sketcherMinimizer::getShapeIndex(vector<sketcherMinimizerPointF> shape,
                                     float floatPosition)
{
    float normalizedF = floatPosition;
    while (normalizedF < 0)
        normalizedF += 1.f;
    while (normalizedF >= 1.f)
        normalizedF -= 1.f;
    int counter = shape.size() * normalizedF;
    return counter;
}

vector<vector<sketcherMinimizerResidue*>>
sketcherMinimizer::groupResiduesInSSEs(
    vector<sketcherMinimizerResidue*> residues)
{
    // divide residues by chain
    map<string, vector<sketcherMinimizerResidue*>> chainsMap;
    for (auto res : residues) {
        chainsMap[res->chain].push_back(res);
    }
    // order each chain by residue number
    for (auto& pair : chainsMap) {
        sort(pair.second.begin(), pair.second.end(),
             [](const sketcherMinimizerResidue* firstRes,
                const sketcherMinimizerResidue* secondRes) {
                 return firstRes->resnum < secondRes->resnum;
             });
    }
    int gap = 3;
    // split chains in smaller chunks whenever more than gap consecutive
    // residues are missing
    vector<vector<sketcherMinimizerResidue*>> returnValue;
    for (auto& pair : chainsMap) {
        vector<sketcherMinimizerResidue*> growingChain;
        for (auto res : pair.second) {
            if (!growingChain.empty() &&
                (res->resnum - growingChain.back()->resnum > gap ||
                 res->chain == " " || res->chain.empty())) {
                returnValue.push_back(growingChain);
                growingChain.clear();
            }
            growingChain.push_back(res);
        }
        if (!growingChain.empty()) {
            returnValue.push_back(growingChain);
        }
    }
    return returnValue;
}

vector<sketcherMinimizerPointF> sketcherMinimizer::shapeAroundLigand(int crownN)
{
    // return crownN-th crown around ligand
    float distanceOfFirstCrown = 60;
    float distanceBetweenCrowns = 60;
    // find limits
    auto atoms = _atoms;
    auto bonds = _bonds;
    float border = distanceBetweenCrowns * crownN + distanceOfFirstCrown;
    float minX = atoms[0]->coordinates.x();
    float maxX = atoms[0]->coordinates.x();
    float minY = atoms[0]->coordinates.y();
    float maxY = atoms[0]->coordinates.y();
    float maxPocketD = 0;
    for (auto atom : atoms) {
        float distance = atom->m_pocketDistance;
        if (distance > maxPocketD) {
            maxPocketD = distance;
        }
        float newX = atom->coordinates.x();
        float newY = atom->coordinates.y();

        if (minX > newX) {
            minX = newX;
        }
        if (maxX < newX) {
            maxX = newX;
        }
        if (minY > newY) {
            minY = newY;
        }
        if (maxY < newY) {
            maxY = newY;
        }
    }
    maxPocketD += 10; // to account for cutoffs at borders
    minX -= border + maxPocketD;
    maxX += border + maxPocketD;
    minY -= border + maxPocketD;
    maxY += border + maxPocketD;

    /*run a marching square algorithm on a grid of x_interval spacing*/
    sketcherMinimizerMarchingSquares ms;
    float x_interval = 20.f;
    ms.initialize(minX, maxX, minY, maxY, x_interval);
    for (unsigned int j = 0; j < ms.getYN(); j++) {
        for (unsigned int i = 0; i < ms.getXN(); i++) {
            float pointX = ms.toRealx(i);
            float pointY = ms.toRealy(j);
            sketcherMinimizerPointF p(pointX, pointY);

            float shortestD = -1;
            for (auto a : atoms) {
                float dist = a->m_pocketDistance;

                auto vect = a->coordinates - p;
                float D = vect.length();
                D -= dist + border;
                if (D < shortestD || shortestD < 0)
                    shortestD = D;
            }

            for (auto b : bonds) {

                sketcherMinimizerPointF sp1(b->startAtom->coordinates);
                sketcherMinimizerPointF sp2(b->endAtom->coordinates);

                float distancePercentage = 1.f;
                float D2 = sketcherMinimizerMaths::squaredDistancePointSegment(
                    p, sp1, sp2, &distancePercentage);
                float D = sqrt(D2);
                float distance2 = b->startAtom->m_pocketDistance;
                float distance1 = b->endAtom->m_pocketDistance;
                float dist = distance1 * distancePercentage +
                             distance2 * (1 - distancePercentage);
                D -= dist + border;

                if (D < shortestD)
                    shortestD = D;
            }

            ms.setValue(shortestD, i, j);
        }
    }
    ms.setThreshold(0);
    ms.run();
    auto result = ms.getOrderedCoordinatesPoints();
    sort(result.begin(), result.end(), [](const vector<float>& firstContour,
                                          const vector<float>& secondContour) {
        return firstContour.size() > secondContour.size();
    });
    vector<sketcherMinimizerPointF> returnValue;
    if (result.size() > 0) {
        for (unsigned int i = 0; i < result.at(0).size(); i += 2) {
            returnValue.push_back(sketcherMinimizerPointF(
                result.at(0).at(i), result.at(0).at(i + 1)));
        }
    }
    return returnValue;
}

float sketcherMinimizer::scoreResiduePosition(
    int index, const vector<sketcherMinimizerPointF>& shape, int shapeN,
    vector<bool>& , sketcherMinimizerResidue* residue)
{
    auto position = shape.at(index);
    float distancePenalty = 0.01f;
    float clashingLigandAtomsPenalty = 100.f;
    vector<sketcherMinimizerAtom*> targets;
    for (auto interactionPartner : residue->residueInteractionPartners) {
        if (interactionPartner->coordinatesSet)
            targets.push_back(interactionPartner);
    }
    float interactionsF = 1.f;
    if (targets.empty()) {
        interactionsF = 0.2f;
        targets.push_back(residue->m_closestLigandAtom);
    }
    float score = 0.f;
    for (auto target : targets) {
        int clashingLigandAtoms = 0;
        for (auto ligandAtom : _atoms) {
            if (ligandAtom == target)
                continue;
            auto ligandAtomPos = ligandAtom->coordinates;
            float squareDist =
                sketcherMinimizerMaths::squaredDistancePointSegment(
                    ligandAtomPos, position, target->coordinates);
            if (squareDist < 40 * 40) {
                clashingLigandAtoms++;
            }
        }
        auto distance = sketcherMinimizerMaths::squaredDistance(
                            target->coordinates, position) -
                        (shapeN * 50) * (shapeN * 50);
        score +=
            interactionsF * (distancePenalty * distance +
                             clashingLigandAtoms * clashingLigandAtomsPenalty);
    }
    return score;
}

void sketcherMinimizer::placeResidues(vector<sketcherMinimizerAtom*> atoms)
{
    if (!_residues.size()) {
        return;
    }
    if (!atoms.size()) {
        placeResiduesProteinOnlyMode();
        return;
    }
    findClosestAtomToResidues(atoms);

    placeResiduesInCrowns();
    m_minimizer.minimizeResidues();
    return;
}

/*
 move mol around to avoid clashes with other already placed molecules. Explore a
 grid of @levels concentric levels, with #gridD resolution. @distanceFromAtoms
 is the minimum clash distance to reject a position.
 */

sketcherMinimizerPointF
sketcherMinimizer::exploreMolPosition(sketcherMinimizerMolecule* mol,
                                      unsigned int levels, float gridD,
                                      float distanceFromAtoms)
{

    sketcherMinimizerPointF v(0, 0), centerOfGrid(0, 0);
    for (unsigned int i = 0; i < levels; i++) {

        vector<sketcherMinimizerPointF> pointstoTest;
        sketcherMinimizerPointF top =
            centerOfGrid + sketcherMinimizerPointF(0.f, (1 + i) * gridD);
        sketcherMinimizerPointF bottom =
            centerOfGrid + sketcherMinimizerPointF(0.f, -((1 + i) * gridD));
        sketcherMinimizerPointF right =
            centerOfGrid + sketcherMinimizerPointF((1 + i) * gridD, 0.f);
        sketcherMinimizerPointF left =
            centerOfGrid + sketcherMinimizerPointF(-((1 + i) * gridD), 0.f);

        pointstoTest.push_back(centerOfGrid);
        pointstoTest.push_back(right);
        pointstoTest.push_back(left);
        pointstoTest.push_back(bottom);
        pointstoTest.push_back(top);

        for (unsigned int j = 0; j < i; j++) {
            pointstoTest.push_back(
                right + sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                right - sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                left + sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                left - sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                bottom + sketcherMinimizerPointF(gridD * (j + 1), 0.f));
            pointstoTest.push_back(
                bottom - sketcherMinimizerPointF(gridD * (j + 1), 0.f));
            pointstoTest.push_back(
                top + sketcherMinimizerPointF(gridD * (j + 1), 0.f));
            pointstoTest.push_back(
                top - sketcherMinimizerPointF(gridD * (j + 1), 0.f));
        }

        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF((1 + i) * gridD, (1 + i) * gridD));
        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF((1 + i) * gridD, -((1 + i) * gridD)));
        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF(-((1 + i) * gridD), (1 + i) * gridD));
        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF(-((1 + i) * gridD), -((1 + i) * gridD)));

        bool noClash = true;
        if (distanceFromAtoms < 0)
            distanceFromAtoms = bondLength * 1.8;
        float dist = distanceFromAtoms;

        for (unsigned int pc = 0; pc < pointstoTest.size(); pc++) {
            noClash = true;
            v = pointstoTest[pc];
            foreach (sketcherMinimizerAtom* at, mol->_atoms) {
                sketcherMinimizerPointF placeNextTo = at->coordinates + v;
                foreach (sketcherMinimizerMolecule* m, _molecules) {
                    if (!m->isPlaced)
                        continue;
                    if (m == mol)
                        continue;
                    foreach (sketcherMinimizerAtom* a, m->_atoms) {
                        dist = distanceFromAtoms;

                        if (((a->coordinates.x() < placeNextTo.x() + dist) &&
                             (a->coordinates.x() > placeNextTo.x() - dist)) &&
                            ((a->coordinates.y() < placeNextTo.y() + dist) &&
                             (a->coordinates.y() > placeNextTo.y() - dist))) {
                            noClash = false;
                            break;
                        }
                    }

                    if (!noClash)
                        break;
                }

                if (!noClash)
                    break;
            }
            if (noClash)
                break;
        }
        if (noClash)
            break;
    }
    return v;
}

sketcherMinimizerPointF sketcherMinimizer::exploreGridAround(
    sketcherMinimizerPointF centerOfGrid, unsigned int levels, float gridD,
    float dx, float dy, float distanceFromAtoms, bool watermap,
    sketcherMinimizerResidue* residueForInteractions,
    sketcherMinimizerPointF direction)
{
    sketcherMinimizerPointF placeNextTo = centerOfGrid;

    for (unsigned int i = 0; i < levels; i++) {

        vector<sketcherMinimizerPointF> pointstoTest;
        sketcherMinimizerPointF top =
            centerOfGrid + sketcherMinimizerPointF(0.f, (1 + i) * gridD);
        sketcherMinimizerPointF bottom =
            centerOfGrid + sketcherMinimizerPointF(0.f, -((1 + i) * gridD));
        sketcherMinimizerPointF right =
            centerOfGrid + sketcherMinimizerPointF((1 + i) * gridD, 0.f);
        sketcherMinimizerPointF left =
            centerOfGrid + sketcherMinimizerPointF(-((1 + i) * gridD), 0.f);

        pointstoTest.push_back(centerOfGrid);
        pointstoTest.push_back(right);
        pointstoTest.push_back(left);
        pointstoTest.push_back(bottom);
        pointstoTest.push_back(top);

        for (unsigned int j = 0; j < i; j++) {
            pointstoTest.push_back(
                right + sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                right - sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                left + sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                left - sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                bottom + sketcherMinimizerPointF(gridD * (j + 1), 0.f));
            pointstoTest.push_back(
                bottom - sketcherMinimizerPointF(gridD * (j + 1), 0.f));
            pointstoTest.push_back(
                top + sketcherMinimizerPointF(gridD * (j + 1), 0.f));
            pointstoTest.push_back(
                top - sketcherMinimizerPointF(gridD * (j + 1), 0.f));
        }

        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF((1 + i) * gridD, (1 + i) * gridD));
        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF((1 + i) * gridD, -((1 + i) * gridD)));
        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF(-((1 + i) * gridD), (1 + i) * gridD));
        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF(-((1 + i) * gridD), -((1 + i) * gridD)));

        bool noClash = true;
        if (distanceFromAtoms < 0)
            distanceFromAtoms = bondLength * 1.8;
        float distanceFromResidues = bondLength * 1.3;
        float watermapDistance = 10;
        float dist = distanceFromAtoms;

        sketcherMinimizerPointF directionNormal(-direction.y(), direction.x());
        directionNormal.normalize();
        for (unsigned int pc = 0; pc < pointstoTest.size(); pc++) {
            noClash = true;
            sketcherMinimizerPointF point = pointstoTest[pc] - centerOfGrid;
            placeNextTo = point.y() * direction + point.x() * directionNormal +
                          centerOfGrid;
            foreach (sketcherMinimizerMolecule* m, _molecules) {
                if (!m->isPlaced)
                    continue;
                foreach (sketcherMinimizerAtom* a, m->_atoms) {
                    if (a->isResidue())
                        dist = distanceFromResidues;
                    else
                        dist = distanceFromAtoms;
                    if (watermap) {
                        if (!a->isResidue())
                            continue;
                        dist = watermapDistance;
                    }
                    if (((a->coordinates.x() < placeNextTo.x() + dist + dx) &&
                         (a->coordinates.x() > placeNextTo.x() - dist - dx)) &&
                        ((a->coordinates.y() < placeNextTo.y() + dist + dy) &&
                         (a->coordinates.y() > placeNextTo.y() - dist - dy))) {
                        noClash = false;
                        break;
                    }
                    if (residueForInteractions) {
                        for (auto partnerOfA : a->residueInteractionPartners) {
                            if (a == residueForInteractions ||
                                partnerOfA == residueForInteractions ||
                                !a->coordinatesSet ||
                                !partnerOfA->coordinatesSet) {
                                continue;
                            }

                            float squareD = sketcherMinimizerMaths::
                                squaredDistancePointSegment(
                                    placeNextTo, a->coordinates,
                                    partnerOfA->coordinates);
                            if (squareD <
                                distanceFromResidues * distanceFromResidues) {
                                noClash = false;
                                break;
                            }

                            for (auto partner :
                                 residueForInteractions
                                     ->residueInteractionPartners) {
                                if (!partner->coordinatesSet) {
                                    continue;
                                }
                                if (sketcherMinimizerMaths::
                                        intersectionOfSegments(
                                            placeNextTo, partner->coordinates,
                                            a->coordinates,
                                            partnerOfA->coordinates)) {
                                    noClash = false;
                                    break;
                                }
                            }
                        }
                    }
                }

                if (!noClash)
                    break;
            }

            if (noClash)
                break;
        }
        if (noClash)
            break;
    }
    return placeNextTo;
}

vector<proximityData> sketcherMinimizer::buildProximityDataVector(
    vector<sketcherMinimizerMolecule*>& proximityMols,
    map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap)
{

    vector<proximityData> proximityDataVector;

    foreach (sketcherMinimizerMolecule* mol, proximityMols) {
        proximityData data;
        sketcherMinimizerAtom* metaAtom = molMap[mol];
        vector<sketcherMinimizerPointF> additionVectors(
            metaAtom->neighbors.size(), sketcherMinimizerPointF(0.f, 0.f));
        vector<sketcherMinimizerPointF> centers(
            metaAtom->neighbors.size(), sketcherMinimizerPointF(0.f, 0.f));
        vector<int> counters(metaAtom->neighbors.size(), 0);
        foreach (sketcherMinimizerBond* pr, mol->m_proximityRelations) {
            sketcherMinimizerAtom* otherMetaAtom = NULL;
            sketcherMinimizerAtom* targetAtom = NULL;
            if (pr->startAtom->molecule == mol &&
                !(pr->endAtom->molecule == mol)) {
                otherMetaAtom = molMap[pr->endAtom->molecule];
                targetAtom = pr->startAtom;
            } else if (pr->endAtom->molecule == mol &&
                       !(pr->startAtom->molecule == mol)) {
                otherMetaAtom = molMap[pr->startAtom->molecule];
                targetAtom = pr->endAtom;
            }
            if (otherMetaAtom) {
                for (unsigned int i = 0; i < metaAtom->neighbors.size(); i++) {
                    if (metaAtom->neighbors[i] == otherMetaAtom) {
                        additionVectors[i] +=
                            targetAtom->getSingleAdditionVector();
                        centers[i] += targetAtom->coordinates;
                        counters[i]++;
                    }
                }
            }
        }
        for (unsigned int i = 0; i < centers.size(); i++) {
            if (counters[i] > 0)
                centers[i] /= counters[i];
            additionVectors[i].normalize();
        }
        data.additionVectors = additionVectors;
        data.centers = centers;
        data.counters = counters;
        proximityDataVector.push_back(data);
    }
    return proximityDataVector;
}

void sketcherMinimizer::rotateMoleculesWithProximityRelations(
    vector<sketcherMinimizerMolecule*>& proximityMols,
    map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap,
    vector<proximityData>& proximityDataVector)
{

    for (unsigned int m = 0; m < proximityMols.size(); m++) {
        sketcherMinimizerMolecule* mol = proximityMols[m];
        sketcherMinimizerAtom* metaAtom = molMap[mol];
        vector<sketcherMinimizerPointF> additionVectors =
            proximityDataVector[m].additionVectors;
        vector<sketcherMinimizerPointF> centers =
            proximityDataVector[m].centers;
        if (mol->_atoms.size() < 2)
            continue;

        sketcherMinimizerPointF direction(1, 0);
        if (metaAtom->bonds.size() == 1) {
            direction =
                metaAtom->coordinates - metaAtom->neighbors[0]->coordinates;
            sketcherMinimizerPointF p1 = additionVectors[0];
            p1 *= -1;
            sketcherMinimizerPointF p3 = direction;
            float rotationAngle = sketcherMinimizerMaths::signedAngle(
                p1, sketcherMinimizerPointF(0, 0), p3);
            rotationAngle *= -M_PI / 180.f;
            float s = sin(rotationAngle);
            float c = cos(rotationAngle);

            foreach (sketcherMinimizerAtom* a, mol->_atoms) {
                sketcherMinimizerPointF coords = a->coordinates - centers[0];
                coords.rotate(s, c);
                a->coordinates = coords + centers[0];
            }
        } else if (metaAtom->bonds.size() > 1) {
            vector<sketcherMinimizerPointF> v1, v2;
            foreach (sketcherMinimizerAtom* n, metaAtom->neighbors) {
                v1.push_back(n->coordinates - metaAtom->coordinates);
            }
            v2 = additionVectors;
            float rotMat[4];
            alignmentMatrix(v1, v2, rotMat);
            sketcherMinimizerPointF center = mol->center();
            foreach (sketcherMinimizerAtom* a, mol->_atoms) {
                sketcherMinimizerPointF coords = a->coordinates - center;
                float x = coords.x();
                float y = coords.y();
                sketcherMinimizerPointF newCoords(x * rotMat[0] + y * rotMat[1],
                                                  x * rotMat[2] +
                                                      y * rotMat[3]);
                a->coordinates = center + newCoords;
            }
        }
    }
}

void sketcherMinimizer::translateMoleculesWithProximityRelations(
    vector<sketcherMinimizerMolecule*>& proximityMols,
    map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap,
    map<sketcherMinimizerMolecule*, sketcherMinimizerPointF>& templateCenters,
    vector<proximityData>& )
{

    // placing
    int counterN = 1;
    bool cleverPlacing = false;
    do {
        cleverPlacing = !cleverPlacing; // alternatively try to be smart aboout
                                        // single atom mols
        if (!cleverPlacing)
            counterN++;

        for (unsigned int m = 0; m < proximityMols.size(); m++) {
            sketcherMinimizerMolecule* mol = proximityMols[m];

            bool residue = false;
            if (mol->_atoms.size() == 1)
                if (mol->_atoms[0]->isResidue())
                    residue = true;
            if (!residue) {
                if (mol->hasConstrainedFragments) {
                    mol->isPlaced = true;
                    continue;
                }
            }
            if (mol->hasFixedFragments) {
                mol->isPlaced = true;
                continue;
            }
            if (mol->m_proximityRelations.size()) {
                sketcherMinimizerPointF atomsCenter =
                    sketcherMinimizerPointF(0, 0);
                int atomsN = 0;
                foreach (sketcherMinimizerBond* pr, mol->m_proximityRelations) {
                    if (pr->startAtom->molecule == mol &&
                        pr->endAtom->molecule != mol) {
                        atomsCenter += pr->startAtom->coordinates;
                        atomsN++;
                    } else if (pr->endAtom->molecule == mol &&
                               pr->startAtom->molecule != mol) {
                        atomsCenter += pr->endAtom->coordinates;
                        atomsN++;
                    }
                }
                if (atomsN > 0)
                    atomsCenter /= atomsN;

                /*positioning */
                sketcherMinimizerPointF placeNextTo = templateCenters[mol];
                placeNextTo *= counterN;

                foreach (sketcherMinimizerAtom* a, mol->_atoms) {
                    a->coordinates += placeNextTo - atomsCenter;
                }

                mol->isPlaced = true;
            }
        }
        if (cleverPlacing) { // replace single terminal atoms
            for (unsigned int m = 0; m < proximityMols.size(); m++) {
                sketcherMinimizerMolecule* mol = proximityMols[m];
                if (mol->_atoms.size() == 1) {
                    sketcherMinimizerAtom* metaAtom = molMap[mol];
                    if (metaAtom->neighbors.size() == 1) {

                        int bondsN = 0;
                        sketcherMinimizerPointF coords(0, 0);
                        foreach (sketcherMinimizerBond* pr,
                                 mol->m_proximityRelations) {
                            if (pr->startAtom->molecule == mol &&
                                pr->endAtom->molecule != mol) {
                                sketcherMinimizerPointF addV =
                                    pr->endAtom->getSingleAdditionVector();
                                if (addV.length() < SKETCHER_EPSILON)
                                    continue;
                                addV.normalize();
                                addV *= bondLength * counterN;
                                coords += pr->endAtom->coordinates + addV;
                                bondsN++;
                            } else if (pr->endAtom->molecule == mol &&
                                       pr->startAtom->molecule != mol) {
                                sketcherMinimizerPointF addV =
                                    pr->startAtom->getSingleAdditionVector();
                                if (addV.length() < SKETCHER_EPSILON)
                                    continue;

                                addV.normalize();
                                addV *= bondLength * counterN;
                                coords += pr->startAtom->coordinates + addV;
                                bondsN++;
                            }
                        }
                        if (bondsN > 0) {
                            coords /= bondsN;

                            mol->_atoms[0]->coordinates = coords;
                        } else { // a suitable addition Vector could not be
                                 // found, try positionings a bondlength to the
                                 // right of a proximity partner
                            foreach (sketcherMinimizerBond* pr,
                                     mol->m_proximityRelations) {
                                if (pr->startAtom->molecule == mol &&
                                    pr->endAtom->molecule != mol) {
                                    mol->_atoms[0]->coordinates =
                                        pr->endAtom->coordinates;
                                    break;
                                } else if (pr->endAtom->molecule == mol &&
                                           pr->startAtom->molecule != mol) {
                                    mol->_atoms[0]->coordinates =
                                        pr->startAtom->coordinates;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

    } while (m_minimizer.findIntermolecularClashes(proximityMols,
                                                   bondLength * 0.5) &&
             counterN < 10);
}

void sketcherMinimizer::placeMoleculesWithProximityRelations(
    vector<sketcherMinimizerMolecule*> proximityMols)
{
    map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*> molMap;

    sketcherMinimizerMolecule* metaMol = new sketcherMinimizerMolecule;
    foreach (sketcherMinimizerMolecule* mol, _molecules) {
        if (mol->m_proximityRelations.size()) {
            sketcherMinimizerAtom* a = new sketcherMinimizerAtom;
            a->molecule = metaMol;
            metaMol->_atoms.push_back(a);
            molMap[mol] = a;
        }
    }

    foreach (sketcherMinimizerBond* b, m_proximityRelations) {
        if (b->startAtom->molecule == b->endAtom->molecule)
            continue;
        sketcherMinimizerAtom* at1 = molMap[b->startAtom->molecule];
        sketcherMinimizerAtom* at2 = molMap[b->endAtom->molecule];
        bool found = false;
        foreach (sketcherMinimizerBond* b, metaMol->_bonds) {
            if ((b->startAtom == at1 && b->endAtom == at2) ||
                (b->startAtom == at2 && b->endAtom == at1)) {
                found = true;
            }
        }
        if (!found) {
            sketcherMinimizerBond* newBond = new sketcherMinimizerBond;
            newBond->startAtom = at1;
            newBond->endAtom = at2;

            metaMol->_bonds.push_back(newBond);
        }
    }
    sketcherMinimizer min(m_minimizer.getPrecision());

    if (metaMol->_atoms.size()) {
        min.m_fragmentBuilder.m_evenAngles = true;
        min.m_minimizer.m_evenAngles = true;

        min.initialize(metaMol);
        min.findFragments();
        min.m_minimizer.buildFromFragments(true);
        min.m_minimizer.avoidClashes();
        min.bestRotation();
        min.maybeFlip();
        min.arrangeMultipleMolecules();
    }

    bool ligandResidueStyle = true; // positions of molecules are determined
                                    // more by a bigger central molecule than by
                                    // a bonding pattern

    for (auto molecule : min._molecules) {
        if (molecule->_rings.size() >
            0) { // if at least three molecules are connected to each other
                 // (i.e. two residues are connected to each other and both to
                 // the ligand) abort the ligandResidue display style)
            ligandResidueStyle = false;
        }
    }
    sketcherMinimizerMolecule* centralMol = proximityMols[0];
    foreach (sketcherMinimizerMolecule* mol, proximityMols) {
        if (mol->m_proximityRelations.size() >
            centralMol->m_proximityRelations.size())
            centralMol = mol;
        else if (mol->m_proximityRelations.size() ==
                     centralMol->m_proximityRelations.size() &&
                 mol->_atoms.size() > centralMol->_atoms.size())
            centralMol = mol;
    }
    if (centralMol->_atoms.size() < MINIMUM_LIGAND_ATOMS)
        ligandResidueStyle = false;
    map<sketcherMinimizerMolecule*, sketcherMinimizerPointF> templateCenters;

    foreach (sketcherMinimizerMolecule* mol, proximityMols) {
        sketcherMinimizerPointF point(0, 0);

        sketcherMinimizerAtom* at = molMap[mol];
        if (at)
            point = at->coordinates;
        templateCenters[mol] = point;
    }

    if (ligandResidueStyle) {
        queue<sketcherMinimizerMolecule*> q;
        map<sketcherMinimizerMolecule*, sketcherMinimizerMolecule*> getParent;

        q.push(centralMol);
        while (q.size()) {
            sketcherMinimizerMolecule* mol = q.front();
            q.pop();
            if (mol->isPlaced)
                continue;
            if (mol == centralMol) {
                mol->isPlaced = true;
            } else {
                sketcherMinimizerMolecule* parent = getParent[mol];
                if (parent != NULL)
                    placeMolResidueLigandStyle(mol, parent);
            }
            foreach (sketcherMinimizerBond* b, mol->m_proximityRelations) {
                if (!b->startAtom->molecule
                         ->isPlaced) { // will place a molecule twice if it has
                                       // two relations with mol. This is safe
                                       // cause of the continue for mol->isPlace
                                       // when looping the second time
                    q.push(b->startAtom->molecule);
                    getParent[b->startAtom->molecule] = mol;
                }
                if (!b->endAtom->molecule->isPlaced) {
                    q.push(b->endAtom->molecule);
                    getParent[b->endAtom->molecule] = mol;
                }
            }
        }

    } else {

        vector<proximityData> proximityDataVector =
            buildProximityDataVector(proximityMols, molMap);

        rotateMoleculesWithProximityRelations(proximityMols, molMap,
                                              proximityDataVector);

        translateMoleculesWithProximityRelations(
            proximityMols, molMap, templateCenters, proximityDataVector);
    }
}

void sketcherMinimizer::placeMolResidueLigandStyle(
    sketcherMinimizerMolecule* mol, sketcherMinimizerMolecule* parent)
{

    int n = 0;
    sketcherMinimizerPointF parentV(0, 0);
    sketcherMinimizerPointF parentAdditionV(0, 0);
    sketcherMinimizerPointF v(0, 0);
    sketcherMinimizerPointF additionV(0,
                                      0); // actually using line to centroid, to
                                          // orient the molecule away from the
                                          // ligand
    sketcherMinimizerPointF cent = mol->center();

    foreach (sketcherMinimizerBond* b, mol->m_proximityRelations) {
        sketcherMinimizerAtom *at = NULL, *parentAt = NULL;
        if (b->startAtom->molecule == parent) {
            parentAt = b->startAtom;
            at = b->endAtom;
        } else if (b->endAtom->molecule == parent) {
            at = b->startAtom;
            parentAt = b->endAtom;
        }
        if (at == NULL || parentAt == NULL)
            continue;
        n++;
        sketcherMinimizerPointF paddV = parentAt->getSingleAdditionVector();
        if (b->isResidueInteraction()) {
            sketcherMinimizerResidueInteraction* ri =
                static_cast<sketcherMinimizerResidueInteraction*>(b);
            if (ri->startAtom->molecule == parent &&
                ri->m_otherStartAtoms.size()) {
                paddV = sketcherMinimizerAtom::getSingleAdditionVector(
                    ri->getAllStartAtoms());
            } else if (ri->endAtom->molecule == parent &&
                       ri->m_otherEndAtoms.size()) {
                paddV = sketcherMinimizerAtom::getSingleAdditionVector(
                    ri->getAllEndAtoms());
            }
        }

        paddV.normalize();
        paddV *= bondLength * 3;
        parentV += parentAt->coordinates;
        parentAdditionV += paddV;
        additionV += at->coordinates - cent;
        v += at->coordinates;
    }
    if (n > 0) {
        v /= n;
        parentV /= n;
        parentAdditionV /= n;
        additionV /= n;
        sketcherMinimizerPointF startingPos = parentV + parentAdditionV;
        startingPos = exploreGridAround(startingPos, 15, 10);

        float angle = sketcherMinimizerMaths::signedAngle(
                          startingPos - parentV, sketcherMinimizerPointF(0, 0),
                          -additionV) /
                      180 * M_PI;
        float s = sin(angle);
        float c = cos(angle);

        foreach (sketcherMinimizerAtom* a, mol->_atoms) {
            a->coordinates -= v;
            a->coordinates.rotate(s, c);
            a->coordinates += startingPos;
            a->coordinates.round();
        }

        flipIfCrossingInteractions(mol);

        sketcherMinimizerPointF avoidClashV = exploreMolPosition(
            mol, 15,
            bondLength *
                0.5); // explore positions on a grid of points to solve clashes
        foreach (sketcherMinimizerAtom* a, mol->_atoms) {
            a->coordinates += avoidClashV;
        }
    }
    mol->isPlaced = true;
}

void sketcherMinimizer::flipIfCrossingInteractions(
    sketcherMinimizerMolecule* mol)
{

    for (unsigned int bb = 0; bb < mol->m_proximityRelations.size() - 1; bb++) {
        bool out = false;
        sketcherMinimizerBond* pr1 = mol->m_proximityRelations[bb];
        if (pr1->startAtom->molecule == pr1->endAtom->molecule) {
            continue;
        }
        if (!(pr1->startAtom->molecule->isPlaced ||
              pr1->startAtom->molecule == mol))
            continue;
        if (!(pr1->endAtom->molecule->isPlaced ||
              pr1->endAtom->molecule == mol))
            continue;

        for (unsigned int bb2 = bb + 1; bb2 < mol->m_proximityRelations.size();
             bb2++) {
            sketcherMinimizerBond* pr2 = mol->m_proximityRelations[bb2];
            if (pr2->startAtom->molecule == pr2->endAtom->molecule) {
                continue;
            }
            if (!(pr2->startAtom->molecule->isPlaced ||
                  pr2->startAtom->molecule == mol))
                continue;
            if (!(pr2->endAtom->molecule->isPlaced ||
                  pr2->endAtom->molecule == mol))
                continue;

            if (sketcherMinimizerMaths::intersectionOfSegments(
                    pr1->startAtom->coordinates, pr1->endAtom->coordinates,
                    pr2->startAtom->coordinates, pr2->endAtom->coordinates)) {
                /* mirror the coordinates */
                sketcherMinimizerAtom* p1 = NULL;
                sketcherMinimizerAtom* p2 = NULL;
                if (pr1->startAtom->molecule == mol)
                    p1 = pr1->startAtom;
                else if (pr1->endAtom->molecule == mol)
                    p1 = pr1->endAtom;

                if (pr2->startAtom->molecule == mol)
                    p2 = pr2->startAtom;
                else if (pr2->endAtom->molecule == mol)
                    p2 = pr2->endAtom;
                if (p1 && p2) {
                    sketcherMinimizerPointF middleP =
                        p1->coordinates + p2->coordinates;
                    middleP *= 0.5;
                    sketcherMinimizerPointF p1p2V =
                        p1->coordinates - p2->coordinates;
                    p1p2V.normalize();

                    foreach (sketcherMinimizerAtom* a, mol->_atoms) {

                        sketcherMinimizerPointF v2 = a->coordinates - middleP;
                        float dot =
                            sketcherMinimizerMaths::dotProduct(p1p2V, v2);
                        sketcherMinimizerPointF parallel = p1p2V;
                        parallel *= dot; // parallel component of v2

                        a->coordinates -= 2 * parallel;
                        a->coordinates.round();
                    }

                    out = true;
                    break;
                }
            }
        }
        if (out)
            break;
    }
}

void sketcherMinimizer::arrangeMultipleMolecules()
{
    for (unsigned int i = 0; i < _residues.size(); i++) { // replace residues
        _residues[i]->coordinatesSet = false;
    }
    if (_molecules.size() > 1) {
        // find centers for molecules bound by proximity relations
        vector<sketcherMinimizerMolecule*> proximityMols;
        foreach (sketcherMinimizerMolecule* mol, _molecules) {
            if (mol->m_proximityRelations.size()) {
                proximityMols.push_back(mol);
            }
        }
        sketcherMinimizerPointF center(0.f, 0.f);
        if (proximityMols.size()) {
            placeMoleculesWithProximityRelations(proximityMols);
        } else {
            int maxI = 0;
            maxI = 0;
            int maxSize = _molecules[0]->_atoms.size();
            for (unsigned int i = 0; i < _molecules.size(); i++) {
                sketcherMinimizerMolecule* m = _molecules[i];
                int size = m->_atoms.size();
                if (size > maxSize) {
                    maxI = i;
                    maxSize = size;
                }
            }

            sketcherMinimizerMolecule* centralMol = _molecules[maxI];
            centralMol->isPlaced = true;
            foreach (sketcherMinimizerAtom* a, _atoms)
                a->_generalUseVisited =
                    false; // using _generalUseVisited to keep track of charged
                           // atoms that have already been used for counterions

            center = centralMol->center();
        }

        // placing non counterions
        bool foundCounterion = true;
        while (foundCounterion) {
            foundCounterion = false;
            foreach (sketcherMinimizerMolecule* mol, _molecules) {
                bool residue = false;
                if (mol->_atoms.size() == 1)
                    if (mol->_atoms[0]->isResidue())
                        residue = true;
                if (!residue) {
                    if (mol->hasConstrainedFragments)
                        mol->isPlaced = true;
                }
                if (mol->hasFixedFragments)
                    mol->isPlaced = true;

                if (mol->isPlaced)
                    continue;
                if (residue)
                    continue;
                int charge = mol->totalCharge();
                if (charge == 0) {
                    sketcherMinimizerPointF counterionMin, counterionMax;
                    sketcherMinimizerPointF placeNextTo = center;
                    mol->boundingBox(counterionMin, counterionMax);
                    float counteriondx =
                        (counterionMax.x() - counterionMin.x()) * .5f;
                    float counteriondy =
                        (counterionMax.y() - counterionMin.y()) * .5f;
                    sketcherMinimizerPointF counterionCenter =
                        (counterionMax + counterionMin) * .5f;

                    foundCounterion = true;
                    // explore a grid around placeNextTo to find a suitable
                    // place

                    sketcherMinimizerPointF centerOfGrid = placeNextTo;

                    float gridD = bondLength;
                    placeNextTo = exploreGridAround(centerOfGrid, 10, gridD,
                                                    counteriondx, counteriondy);

                    foreach (sketcherMinimizerAtom* a, mol->_atoms) {
                        a->coordinates += placeNextTo - counterionCenter;
                    }

                    mol->isPlaced = true;
                }
            }
        }

        // placing counterions
        foundCounterion = true;
        while (foundCounterion) {
            foundCounterion = false;
            foreach (sketcherMinimizerMolecule* mol, _molecules) {
                if (mol->isPlaced)
                    continue;
                int charge = mol->totalCharge();
                if (charge != 0) {
                    sketcherMinimizerPointF counterionMin, counterionMax;
                    sketcherMinimizerPointF placeNextTo = center;
                    mol->boundingBox(counterionMin, counterionMax);
                    float counteriondx =
                        (counterionMax.x() - counterionMin.x()) * .5f;
                    float counteriondy =
                        (counterionMax.y() - counterionMin.y()) * .5f;
                    sketcherMinimizerPointF counterionCenter =
                        (counterionMax + counterionMin) * .5f;

                    foundCounterion = true;

                    // find an already placed charged atom to place the
                    // counterion next to

                    foreach (sketcherMinimizerMolecule* m, _molecules) {
                        bool found = false;
                        if (!m->isPlaced)
                            continue;
                        foreach (sketcherMinimizerAtom* a, m->_atoms) {
                            if (a->charge == 0)
                                continue;
                            if (a->_generalUseVisited)
                                continue;
                            if (a->charge * charge < 0) {
                                a->_generalUseVisited = true;
                                placeNextTo = a->coordinates;
                                found = true;
                                break;
                            }
                        }
                        if (found)
                            break;
                    }

                    // explore a grid around placeNextTo to find a suitable
                    // place

                    sketcherMinimizerPointF centerOfGrid = placeNextTo;
                    float gridD = bondLength;

                    placeNextTo =
                        exploreGridAround(centerOfGrid, 10, gridD, counteriondx,
                                          counteriondy, bondLength * 0.8);

                    foreach (sketcherMinimizerAtom* a, mol->_atoms) {
                        a->coordinates += placeNextTo - counterionCenter;
                    }

                    mol->isPlaced = true;
                }
            }
        }
        vector<sketcherMinimizerAtom*> ligandAtoms;
        foreach (sketcherMinimizerAtom* a, _atoms)
            if (a->m_isLigand)
                ligandAtoms.push_back(a);
        placeResidues(ligandAtoms);
    }
}

void sketcherMinimizer::initializeFragments()
{

    if (!_fragments.size()) {
        cerr << "Sketcherlibs warning: no fragments to initialize" << endl;
        return;
    }

    foreach (sketcherMinimizerFragment* indf, _independentFragments) {
        assignNumberOfChildrenAtomsFromHere(
            indf); // recursively assign it to children
    }
    foreach (sketcherMinimizerFragment* f, _fragments) {
        m_fragmentBuilder.initializeCoordinates(f);
    }

    foreach (sketcherMinimizerFragment* indf, _independentFragments) {
        assignLongestChainFromHere(indf); // recursively assign it to children
    }
}

bool sketcherMinimizer::alignWithParentDirectionConstrained(
    sketcherMinimizerFragment* fragment, sketcherMinimizerPointF position,
    float angle)
{
    vector<sketcherMinimizerPointF> templates, plainCoordinates,
        flippedCoordinates;
    float sine = sin(angle);
    float cosine = cos(angle);
    for (auto atom : fragment->_coordinates) {
        if (atom.first->constrained) {
            sketcherMinimizerPointF plainCoordinatesAtom = atom.second;
            sketcherMinimizerPointF flippedCoordinatesAtom(
                plainCoordinatesAtom.x(), -plainCoordinatesAtom.y());
            plainCoordinatesAtom.rotate(sine, cosine);
            flippedCoordinatesAtom.rotate(sine, cosine);
            templates.push_back(atom.first->templateCoordinates);
            plainCoordinates.push_back(plainCoordinatesAtom + position);
            flippedCoordinates.push_back(flippedCoordinatesAtom + position);
        }
    }
    float scorePlain =
        roundToTwoDecimalDigits(RMSD(templates, plainCoordinates));
    float scoreFlipped =
        roundToTwoDecimalDigits(RMSD(templates, flippedCoordinates));
    return (scoreFlipped < scorePlain);
}

vector<sketcherMinimizerBond*>
sketcherMinimizer::getAllTerminalBonds(sketcherMinimizerFragment* fragment)
{
    vector<sketcherMinimizerBond*> bonds;
    for (auto bond : fragment->getBonds()) {
        if (bond->isResidueInteraction())
            continue;
        if (bond->startAtom->neighbors.size() == 1 ||
            bond->endAtom->neighbors.size() == 1) {
            bonds.push_back(bond);
        }
    }
    for (auto child : fragment->_children) {
        bonds.push_back(child->_bondToParent);
    }
    if (fragment->getParent()) {
        bonds.push_back(fragment->_bondToParent);
    }
    return bonds;
}

vector<std::pair<sketcherMinimizerPointF, float>>
sketcherMinimizer::findDirectionsToAlignWith(
    sketcherMinimizerFragment* fragment)
{
    vector<std::pair<sketcherMinimizerPointF, float>> chainDirs;

    sketcherMinimizerPointF origin =
        (fragment->_bondToParent->startAtom->coordinates +
         fragment->_bondToParent->endAtom->coordinates) *
        0.5;
    vector<sketcherMinimizerBond*> parentEndBonds =
        getAllTerminalBonds(fragment->getParent());
    for (auto bond : parentEndBonds) {
        if (bond->endAtom->fragment == fragment)
            continue;
        sketcherMinimizerPointF direction =
            origin -
            (bond->startAtom->coordinates + bond->endAtom->coordinates) * 0.5;
        direction.normalize();
        float score = 1.f;
        if (bond->bondOrder == 2)
            score *= SCORE_MULTIPLIER_FOR_DOUBLE_BONDS;
        if ((bond->startAtom->neighbors.size() == 1 &&
             bond->startAtom->atomicNumber != 6) ||
            (bond->endAtom->neighbors.size() == 1 &&
             bond->endAtom->atomicNumber != 6)) {
            score *= SCORE_MULTIPLIER_FOR_SINGLE_BONDED_HETEROATOMS;
        }
        if (bond->endAtom->fragment != fragment->getParent() ||
            bond->startAtom->fragment != fragment->getParent()) {
            score = bond->endAtom->fragment->longestChainFromHere *
                    SCORE_MULTIPLIER_FOR_FRAGMENTS;
            if (fragment->getParent()->getParent() &&
                bond->startAtom->fragment ==
                    fragment->getParent()->getParent()) {
                score *= 100;
            }
        }
        chainDirs.push_back(
            std::pair<sketcherMinimizerPointF, float>(direction, score));
    }
    return chainDirs;
}

float sketcherMinimizer::testAlignment(
    sketcherMinimizerPointF direction,
    std::pair<sketcherMinimizerPointF, float> templat)
{
    float dot = sketcherMinimizerMaths::dotProduct(direction, templat.first);
    if (dot < 0)
        dot = 0;
    float score = dot * dot;
    if (dot > 1 - SKETCHER_EPSILON)
        score += 1000;
    score *= templat.second;
    return score;
}

sketcherMinimizerPointF sketcherMinimizer::scoreDirections(
    sketcherMinimizerFragment* fragment, float angle,
    vector<std::pair<sketcherMinimizerPointF, float>> directions, bool& invert)
{
    float sine = sin(angle);
    float cosine = cos(angle);
    float bestScore = 0.f;
    sketcherMinimizerPointF bestDirection(1.f, 0.f);
    vector<sketcherMinimizerBond*> terminalBonds =
        getAllTerminalBonds(fragment);
    for (auto bond : terminalBonds) {
        if (bond->startAtom->fragment != fragment)
            continue;
        sketcherMinimizerPointF bondDirectionPlain =
            (fragment->_coordinates[bond->startAtom] +
             fragment->_coordinates[bond->endAtom]) *
                0.5 -
            sketcherMinimizerPointF(-bondLength * 0.5, 0);
        bondDirectionPlain.normalize();
        sketcherMinimizerPointF bondDirectionInverted(bondDirectionPlain.x(),
                                                      -bondDirectionPlain.y());
        bondDirectionPlain.rotate(sine, cosine);
        bondDirectionInverted.rotate(sine, cosine);
        float scoreModifier = 1.f;
        if (bond->bondOrder == 2)
            scoreModifier *= SCORE_MULTIPLIER_FOR_DOUBLE_BONDS;
        if ((bond->startAtom->neighbors.size() == 1 &&
             bond->startAtom->atomicNumber != 6) ||
            (bond->endAtom->neighbors.size() == 1 &&
             bond->endAtom->atomicNumber != 6)) {
            scoreModifier *= SCORE_MULTIPLIER_FOR_SINGLE_BONDED_HETEROATOMS;
        }
        if (bond->endAtom->fragment != fragment) {
            scoreModifier = bond->endAtom->fragment->longestChainFromHere *
                            SCORE_MULTIPLIER_FOR_FRAGMENTS;
        }
        for (auto direction : directions) {
            float scorePlain =
                testAlignment(bondDirectionPlain, direction) * scoreModifier;
            if (scorePlain > bestScore) {
                bestScore = scorePlain;
                bestDirection = direction.first;
                invert = false;
            }
            float scoreInverted =
                testAlignment(bondDirectionInverted, direction) * scoreModifier;
            if (scoreInverted > bestScore) {
                bestScore = scoreInverted;
                bestDirection = direction.first;
                invert = true;
            }
        }
    }
    return bestDirection;
}

bool sketcherMinimizer::alignWithParentDirectionUnconstrained(
    sketcherMinimizerFragment* fragment, float angle)
{
    vector<std::pair<sketcherMinimizerPointF, float>> directions =
        findDirectionsToAlignWith(fragment);
    bool invert = false;
    scoreDirections(fragment, angle, directions, invert);
    return invert;
}

void sketcherMinimizer::alignWithParentDirection(
    sketcherMinimizerFragment* f, sketcherMinimizerPointF position, float angle)
{
    // deciding which "side" the fragment will be drawn, rotating 180° around
    // the axis of its bond to parent
    if (f->fixed)
        return;
    bool invert = (f->constrained
                       ? alignWithParentDirectionConstrained(f, position, angle)
                       : alignWithParentDirectionUnconstrained(f, angle));
    if (invert) {
        for (auto& atom : f->_coordinates) {
            atom.second.setY(-atom.second.y());
        }
        for (auto atom : f->getAtoms()) {
            if (atom->hasStereochemistrySet) {
                for (auto bond : atom->bonds) {
                    bond->isWedge = !bond->isWedge;
                }
            }
        }
    }
}

void sketcherMinimizer::assignNumberOfChildrenAtomsFromHere(
    sketcherMinimizerFragment* f)
{
    float cumulatedNumberOfAtoms = 0;
    float cumulatedNumberOfAtomsRanks = 0;
    float childrenAtoms = 0;
    foreach (sketcherMinimizerFragment* child, f->_children) {
        assignNumberOfChildrenAtomsFromHere(child);
        cumulatedNumberOfAtoms += child->numberOfChildrenAtoms;
        cumulatedNumberOfAtomsRanks += child->numberOfChildrenAtomsRank;
        childrenAtoms += child->getAtoms().size();
    }
    f->numberOfChildrenAtoms = cumulatedNumberOfAtoms + childrenAtoms;
    f->numberOfChildrenAtomsRank =
        0.01f * cumulatedNumberOfAtomsRanks + childrenAtoms;
}

void sketcherMinimizer::assignLongestChainFromHere(sketcherMinimizerFragment* f)
{
    float longestDist = 0;
    foreach (sketcherMinimizerFragment* child, f->_children) {
        assignLongestChainFromHere(child);
        if (child->longestChainFromHere > longestDist) {
            longestDist = child->longestChainFromHere;
        }
    }
    sketcherMinimizerPointF positionFromParent(0.f, 0.f);
    if (f->getParent()) {
        positionFromParent =
            f->getParent()->_coordinates[f->_bondToParent->endAtom];
    }
    f->longestChainFromHere = longestDist + positionFromParent.length();
}

sketcherMinimizerBond*
sketcherMinimizer::getBond(const sketcherMinimizerAtom* a1,
                           const sketcherMinimizerAtom* a2)
{
    for (unsigned int i = 0; i < a1->neighbors.size(); i++) {
        if (a1->neighbors[i] == a2)
            return a1->bonds[i];
    }
    return NULL;
}

sketcherMinimizerAtom*
sketcherMinimizer::pickBestAtom(vector<sketcherMinimizerAtom*>& atoms)
{

    vector<sketcherMinimizerAtom *> candidates, oldCandidates;

    int biggestSize = atoms[0]->fragment->numberOfChildrenAtoms;
    foreach (sketcherMinimizerAtom* a, atoms) {
        int size = a->fragment->numberOfChildrenAtoms;
        if (size == biggestSize) {
            candidates.push_back(a);
        } else if (size > biggestSize) {
            biggestSize = size;
            candidates.clear();
            candidates.push_back(a);
        }
    }
    if (candidates.size() == 1)
        return candidates[0];
    oldCandidates = candidates;
    candidates.clear();

    biggestSize = oldCandidates[0]->fragment->numberOfChildrenAtomsRank;
    foreach (sketcherMinimizerAtom* a, oldCandidates) {
        int size = a->fragment->numberOfChildrenAtomsRank;
        if (size == biggestSize) {
            candidates.push_back(a);
        } else if (size > biggestSize) {
            biggestSize = size;
            candidates.clear();
            candidates.push_back(a);
        }
    }
    if (candidates.size() == 1)
        return candidates[0];
    oldCandidates = candidates;
    candidates.clear();

    biggestSize = oldCandidates[0]->atomicNumber;
    foreach (sketcherMinimizerAtom* a, oldCandidates) {
        int size = a->atomicNumber;
        if (size == biggestSize) {
            candidates.push_back(a);
        } else if (size > biggestSize) {
            biggestSize = size;
            candidates.clear();
            candidates.push_back(a);
        }
    }
    if (candidates.size() == 1)
        return candidates[0];
    oldCandidates = candidates;
    candidates.clear();

    // give up
    return oldCandidates[0];
}

void sketcherMinimizer::constrainAllAtoms()
{
    //   cerr << "sketcherMinimizer::constrainAllAtoms ()"<<endl;
    foreach (sketcherMinimizerAtom* a, _atoms)
        a->constrained = true;
}

void sketcherMinimizer::constrainAtoms(vector<bool> constrained)
{
    if (constrained.size() == _atoms.size()) {
        for (unsigned int i = 0; i < constrained.size(); i++) {
            if (constrained[i])
                _atoms[i]->constrained = true;
        }
    } else {
        cerr << "warning, wrong size of vector for constrained atoms. Ignoring"
             << endl;
    }
}

void sketcherMinimizer::fixAtoms(vector<bool> fixed)
{
    if (fixed.size() == _atoms.size()) {
        for (unsigned int i = 0; i < fixed.size(); i++) {
            if (fixed[i])
                _atoms[i]->fixed = true;
        }
    } else {
        cerr << "warning, wrong size of vector for fixed atoms. Ignoring"
             << endl;
    }
}

void sketcherMinimizer::findClosestAtomToResidues(
    vector<sketcherMinimizerAtom*> atoms)
{
    if (!atoms.size())
        atoms = _atoms;
    foreach (sketcherMinimizerAtom* r, _residues) {
        float squareD = 9999999.f;
        sketcherMinimizerAtom* closestA = NULL;
        foreach (sketcherMinimizerAtom* a, atoms) {
            if (!a->isResidue()) {
                float diffx = a->m_x3D - r->m_x3D;
                float diffy = a->m_y3D - r->m_y3D;
                float diffz = a->m_z3D - r->m_z3D;

                float newSquareD =
                    diffx * diffx + diffy * diffy + diffz * diffz;
                if (newSquareD < squareD) {
                    squareD = newSquareD;
                    closestA = a;
                }
            }
        }
        static_cast<sketcherMinimizerResidue*>(r)->m_closestLigandAtom =
            closestA;
        if (!r->m_isClashing)
            r->m_isClashing = (squareD < RESIDUE_CLASH_DISTANCE_SQUARED);
    }
    foreach (sketcherMinimizerBond* b, _bonds) {
        if (b->startAtom->isResidue()) {
            static_cast<sketcherMinimizerResidue*>(b->startAtom)
                ->m_closestLigandAtom = b->endAtom;
        }
        if (b->endAtom->isResidue()) {
            static_cast<sketcherMinimizerResidue*>(b->endAtom)
                ->m_closestLigandAtom = b->startAtom;
        }
    }
}

float sketcherMinimizer::RMSD(vector<sketcherMinimizerPointF> templates,
                              vector<sketcherMinimizerPointF> points)
{
    assert(templates.size() == points.size());
    int counter = templates.size();
    float total = 0.f;
    for (unsigned int i = 0; i < templates.size(); i++) {
        //        cerr << templates[i].x () << ", "<< templates[i].y () << "
        //        " << points[i].x()<<", "<<points[i].y ()<<endl;
        sketcherMinimizerPointF diff = templates[i] - points[i];
        total += diff.x() * diff.x() + diff.y() * diff.y();
    }
    if (counter > 0)
        total /= counter;
    return sqrt(total);
}

void sketcherMinimizer::alignmentMatrix(vector<sketcherMinimizerPointF> ref,
                                        vector<sketcherMinimizerPointF> points,
                                        float* m)
{
    float U[4];
    float Sig[4];

    float V[4];
    float a[4];

    a[0] = 0.f;
    a[1] = 0.f;
    a[2] = 0.f;
    a[3] = 0.f;
    assert(ref.size() == points.size());
    for (unsigned int i = 0; i < ref.size(); i++) {
        a[0] += ref[i].x() * points[i].x();
        a[1] += ref[i].y() * points[i].x();
        a[2] += ref[i].x() * points[i].y();
        a[3] += ref[i].y() * points[i].y();
    }
    svd(a, U, Sig, V);
    m[0] = V[0] * U[0] + V[1] * U[1];
    m[1] = V[0] * U[2] + V[1] * U[3];
    m[2] = V[2] * U[0] + V[3] * U[1];
    m[3] = V[2] * U[2] + V[3] * U[3];
}

void sketcherMinimizer::svd(float* a, float* U, float* Sig, float* V)
{
    float a1[4];
    a1[0] = a[0];
    a1[1] = a[2];
    a1[2] = a[1];
    a1[3] = a[3];

    float Su[4];
    Su[0] = a[0] * a1[0] + a[1] * a1[2];
    Su[1] = a[0] * a1[1] + a[1] * a1[3];
    Su[2] = a[2] * a1[0] + a[3] * a1[2];
    Su[3] = a[2] * a1[1] + a[3] * a1[3];

    float phi = 0.5 * atan2(Su[1] + Su[2], Su[0] - Su[3]);
    float cphi = cos(phi);
    cphi = roundToTwoDecimalDigits(cphi);
    float sphi = sin(phi);
    sphi = roundToTwoDecimalDigits(sphi);

    U[0] = cphi * (-1);
    U[1] = -sphi;
    U[2] = sphi * (-1);
    U[3] = cphi;

    float Sw[4];
    Sw[0] = a1[0] * a[0] + a1[1] * a[2];
    Sw[1] = a1[0] * a[1] + a1[1] * a[3];
    Sw[2] = a1[2] * a[0] + a1[3] * a[2];
    Sw[3] = a1[2] * a[1] + a1[3] * a[3];

    float theta = 0.5 * atan2(Sw[1] + Sw[2], Sw[0] - Sw[3]);
    float ctheta = cos(theta);
    float stheta = sin(theta);

    float W[4];
    W[0] = ctheta;
    W[1] = -stheta;
    W[2] = stheta;
    W[3] = ctheta;

    float SUsum = Su[0] + Su[3];
    float SUdif = sqrt((Su[0] - Su[3]) * (Su[0] - Su[3]) + 4 * Su[1] * Su[2]);

    Sig[0] = sqrt((SUsum + SUdif) * 0.5);
    Sig[1] = 0.f;
    Sig[2] = 0.f;
    Sig[3] = sqrt((SUsum - SUdif) * 0.5);

    float U1[4];
    U1[0] = U[0];
    U1[1] = U[2];
    U1[2] = U[1];
    U1[3] = U[3];

    float U1A[4];

    U1A[0] = U1[0] * a[0] + U1[1] * a[2];
    U1A[1] = U1[0] * a[1] + U1[1] * a[3];
    U1A[2] = U1[2] * a[0] + U1[3] * a[2];
    U1A[3] = U1[2] * a[1] + U1[3] * a[3];

    float S[4];

    S[0] = U1A[0] * W[0] + U1A[1] * W[2];
    S[1] = U1A[0] * W[1] + U1A[1] * W[3];
    S[2] = U1A[2] * W[0] + U1A[3] * W[2];
    S[3] = U1A[2] * W[1] + U1A[3] * W[3];
    S[0] = roundToTwoDecimalDigits(S[0]);
    S[1] = roundToTwoDecimalDigits(S[1]);
    S[2] = roundToTwoDecimalDigits(S[2]);
    S[3] = roundToTwoDecimalDigits(S[3]);

    float sign1 = 1.f, sign2 = 1.f;
    if (S[0] < 0.f)
        sign1 = -1.f;
    if (S[3] < 0.f)
        sign2 = -1.f;
    float C[4];

    C[0] = sign1;
    C[1] = 0.f;
    C[2] = 0.f;
    C[3] = sign2;

    V[0] = W[0] * C[0] + W[1] * C[2];
    V[1] = W[0] * C[1] + W[1] * C[3];
    V[2] = W[2] * C[0] + W[3] * C[2];
    V[3] = W[2] * C[1] + W[3] * C[3];
    V[0] = roundToTwoDecimalDigits(V[0]);
    V[1] = roundToTwoDecimalDigits(V[1]);
    V[2] = roundToTwoDecimalDigits(V[2]);
    V[3] = roundToTwoDecimalDigits(V[3]);
}

bool sketcherMinimizer::compare(vector<sketcherMinimizerAtom*> atoms,
                                vector<sketcherMinimizerBond*> bonds,
                                sketcherMinimizerMolecule* templ,
                                vector<unsigned int>& mapping)
{
    if (atoms.size() != templ->_atoms.size())
        return false;
    vector<int> molScores, tempScores;
    int molIter = morganScores(atoms, bonds, molScores);
    int templateIter = morganScores(templ->_atoms, templ->_bonds, tempScores);

    if (molIter != templateIter)
        return false;

    unsigned int size = atoms.size();
    vector<bool> matrix(size * size, false);

    vector<sketcherMinimizerPointF> templateCoordinates;
    vector<vector<int>> molBonds;
    vector<vector<int>> templateBonds;

    //   vector < vector < int > > templateAllowedZChains; //for double bond
    //   chirality
    vector<vector<int>> molCisTransChains;
    vector<bool> molIsCis;

    for (unsigned int ma = 0; ma < size; ma++) {
        vector<int> vec;
        molBonds.push_back(vec);
    }

    for (unsigned int ta = 0; ta < size; ta++) {
        vector<int> vec;
        templateBonds.push_back(vec);
    }

    for (unsigned int ta = 0; ta < size; ta++) {
        templateCoordinates.push_back(templ->_atoms[ta]->coordinates);
    }
    for (unsigned int mb = 0; mb < bonds.size(); mb++) {
        if (bonds[mb]->bondOrder == 2) {

            sketcherMinimizerBond* b = bonds[mb];
            if (b->m_ignoreZE)
                continue;

            bool ok = false;
            if (b->rings.size() == 0)
                ok = true;
            else {
                ok = true;
                for (unsigned int rr = 0; rr < b->rings.size(); rr++) {
                    if (b->rings[rr]->_atoms.size() < MACROCYCLE) {
                        ok = false;
                    }
                }
            }
            if (!ok)
                continue;
            sketcherMinimizerAtom* startA = b->startAtomCIPFirstNeighbor();
            sketcherMinimizerAtom* endA = b->endAtomCIPFirstNeighbor();

            sketcherMinimizerAtom* startN = NULL;
            sketcherMinimizerAtom* endN = NULL;
            for (unsigned int i = 0; i < b->startAtom->neighbors.size(); i++) {
                if (b->startAtom->neighbors[i] == b->endAtom)
                    continue;
                bool found = false;
                for (unsigned int aa = 0; aa < atoms.size(); aa++) {
                    if (atoms[aa] == b->startAtom->neighbors[i]) {
                        startN = b->startAtom->neighbors[i];
                        found = true;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }
            for (unsigned int i = 0; i < b->endAtom->neighbors.size(); i++) {
                if (b->endAtom->neighbors[i] == b->startAtom)
                    continue;
                bool found = false;
                for (unsigned int aa = 0; aa < atoms.size(); aa++) {
                    if (atoms[aa] == b->endAtom->neighbors[i]) {
                        endN = b->endAtom->neighbors[i];
                        found = true;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }

            if (startN && endN && startA && endA) {
                bool isCis = b->isZ;
                if (startN != startA)
                    isCis = !isCis;
                if (endN != endA)
                    isCis = !isCis;
                vector<int> chain;
                chain.push_back(startN->_generalUseN);
                chain.push_back(b->startAtom->_generalUseN);
                chain.push_back(b->endAtom->_generalUseN);
                chain.push_back(endN->_generalUseN);
                molCisTransChains.push_back(chain);
                molIsCis.push_back(isCis);
                // cerr << "chiral double bond" <<isCis<<"    "<<b->isZ<<endl;
            }
        }
    }
    // assuming that _generalUseN is set as the index of each atom
    for (unsigned int mb = 0; mb < bonds.size(); mb++) {
        int in1 = bonds[mb]->startAtom->_generalUseN;
        int in2 = bonds[mb]->endAtom->_generalUseN;

        if (in1 < in2) {
            molBonds[in2].push_back(in1);
        } else {
            molBonds[in1].push_back(in2);
        }
    }
    for (unsigned int tb = 0; tb < templ->_bonds.size(); tb++) {
        int in1 = templ->_bonds[tb]->startAtom->_generalUseN;
        int in2 = templ->_bonds[tb]->endAtom->_generalUseN;
        if (in1 < in2) {
            templateBonds[in2].push_back(in1);
        } else {
            templateBonds[in1].push_back(in2);
        }
    }
    for (unsigned int ma = 0; ma < atoms.size(); ma++) {
        for (unsigned int ta = 0; ta < templ->_atoms.size(); ta++) {
            if (molScores[ma] == tempScores[ta]) {
                matrix[ma * size + ta] = true;
            }
        }
    }
    bool found = false;
    vector<unsigned int> solution;
    for (unsigned int i = 0; i < size; i++) {
        if (!matrix[i])
            continue;
        checkIdentity(solution, i, matrix, templateCoordinates, molBonds,
                      templateBonds, molCisTransChains, molIsCis, size, found,
                      mapping);
        if (found)
            break;
    }
    return found;
}

void sketcherMinimizer::checkIdentity(
    vector<unsigned int> solution, int newSol, vector<bool>& matrix,
    vector<sketcherMinimizerPointF>& templateCoordinates,
    vector<vector<int>>& molBonds, vector<vector<int>>& templateBonds,
    vector<vector<int>>& molCisTransChains, vector<bool>& molIsCis,
    unsigned int size, bool& found, vector<unsigned int>& mapping)
{
    solution.push_back(newSol);
    if (solution.size() == size) {
        // check double bonds stereo
        bool doubleBondsStereoOk = true;

        for (unsigned int i = 0; i < molCisTransChains.size(); i++) {
            sketcherMinimizerPointF p1 =
                templateCoordinates[solution[molCisTransChains[i][0]]];
            sketcherMinimizerPointF p2 =
                templateCoordinates[solution[molCisTransChains[i][1]]];
            sketcherMinimizerPointF p3 =
                templateCoordinates[solution[molCisTransChains[i][2]]];
            sketcherMinimizerPointF p4 =
                templateCoordinates[solution[molCisTransChains[i][3]]];

            if (molIsCis[i] !=
                sketcherMinimizerMaths::sameSide(p1, p4, p2, p3)) {
                doubleBondsStereoOk = false;
                break;
            }
        }

        if (doubleBondsStereoOk) {
            found = true;
            mapping = solution;
        }
    } else {
        for (unsigned int i = 0; i < size; i++) {
            if (found)
                break;
            if (!(matrix[solution.size() * size + i]))
                continue;
            bool check = true;
            for (unsigned int ss = 0; ss < solution.size(); ss++) {
                if (solution[ss] == i) {
                    check = false;
                    break;
                }
            }
            if (check) {
                for (unsigned int bi = 0; bi < molBonds[solution.size()].size();
                     bi++) {
                    check = false;
                    int high = i;
                    int low = solution[molBonds[solution.size()][bi]];
                    if (low > high) {
                        int swap = low;
                        low = high;
                        high = swap;
                    }
                    for (unsigned int ch = 0; ch < templateBonds[high].size();
                         ch++) {
                        if (templateBonds[high][ch] == low) {
                            check = true;
                            break;
                        }
                    }
                    if (check == false)
                        break;
                }
            }

            if (check)
                checkIdentity(solution, i, matrix, templateCoordinates,
                              molBonds, templateBonds, molCisTransChains,
                              molIsCis, size, found, mapping);
        }
    }
}

void sketcherMinimizer::setTemplateFileDir(string dir)
{
    sketcherMinimizer::m_templates.setTemplateDir(dir);
}


static string getTempFileProjDir()
{
    return sketcherMinimizer::m_templates.getTemplateDir();
}

static string getUserTemplateFileName()
{
    const string suffix = "user_templates.mae"; // whatever you wish
    return getTempFileProjDir() + suffix;
}

static void loadTemplate(const string& filename,
                         vector<sketcherMinimizerMolecule*>& templates)
{
    std::ifstream ss(filename);
    schrodinger::mae::Reader r(ss);


    std::shared_ptr<schrodinger::mae::Block> b;
    while ((b = r.next("f_m_ct")) != nullptr) {
        auto molecule = new sketcherMinimizerMolecule();
        int atomCounter = 0;
        // Atom data is in the m_atom indexed block
        {
            const auto atom_data = b->getIndexedBlock("m_atom");
            // All atoms are gauranteed to have these three field names:
            const auto atomic_numbers = atom_data->getIntProperty("i_m_atomic_number");
            const auto xs = atom_data->getRealProperty("r_m_x_coord");
            const auto ys = atom_data->getRealProperty("r_m_y_coord");
            const auto size = atomic_numbers->size();

            // atomic numbers, and x, y, and z coordinates
            for (size_t i=0; i<size; ++i) {
                auto atom = new sketcherMinimizerAtom();
                atom->coordinates = sketcherMinimizerPointF (xs->at(i), ys->at(i));
                atom->atomicNumber = atomic_numbers->at(i);
                atom->_generalUseN = atomCounter++;
                molecule->_atoms.push_back(atom);
            }
        }


        // Bond data is in the m_bond indexed block
        {
            const auto bond_data = b->getIndexedBlock("m_bond");
            // All bonds are gauranteed to have these three field names:
            auto from_atoms = bond_data->getIntProperty("i_m_from");
            auto to_atoms = bond_data->getIntProperty("i_m_to");
            auto orders = bond_data->getIntProperty("i_m_order");
            const auto size = from_atoms->size();

            for (size_t i=0; i<size; ++i) {
                // Maestro atoms are 1 indexed!
                const auto from_atom = from_atoms->at(i) - 1;
                const auto to_atom = to_atoms->at(i) - 1;
                const auto order = orders->at(i);

                auto bond = new sketcherMinimizerBond();
                bond->startAtom = molecule->_atoms.at(from_atom);
                bond->endAtom = molecule->_atoms.at(to_atom);
                bond->bondOrder = order;
                molecule->_bonds.push_back(bond);

            }
        }

        templates.push_back(molecule);
    }

    for (auto mol : templates) {
        // normalize bond length
        vector<float> dds;
        vector<int> ns;
        for (unsigned int i = 0; i < mol->_bonds.size(); i++) {
            sketcherMinimizerPointF v = mol->_bonds[i]->startAtom->coordinates -
            mol->_bonds[i]->endAtom->coordinates;
            float dd = v.x() * v.x() + v.y() * v.y();
            bool found = false;
            for (unsigned int j = 0; j < dds.size(); j++) {
                if (dd * 0.9 < dds[j] && dd * 1.1 > dds[j]) {
                    ns[j]++;
                    found = true;
                    break;
                }
            }
            if (!found) {
                dds.push_back(dd);
                ns.push_back(1);
            }
        }

        if (dds.size()) {
            int maxI = 0;
            for (unsigned int i = 0; i < ns.size(); i++) {
                if (ns[i] > ns[maxI]) {
                    maxI = i;
                }
            }

            float f = sqrt(dds[maxI]);
            for (unsigned int i = 0; i < mol->_atoms.size(); i++) {
                mol->_atoms[i]->coordinates /= f;
            }
        }
    }





}

void sketcherMinimizer::loadTemplates()
{
    static int loaded = 0;
    if (loaded || m_templates.getTemplates().size())
        return;
    string filename =getTempFileProjDir() + "templates.mae";
    loadTemplate(filename, m_templates.getTemplates());

    filename = getUserTemplateFileName();
   // if (ChmFileExists(filename.c_str()))
        loadTemplate(filename, m_templates.getTemplates());

    loaded = 1;
}

int sketcherMinimizer::morganScores(vector<sketcherMinimizerAtom*> atoms,
                                    vector<sketcherMinimizerBond*> bonds,
                                    vector<int>& oldScores)
{
    // assumes that atom[i]->_generalUseN = i
    if (atoms.size() < 2)
        return 0;
    oldScores = vector<int>(atoms.size(), 1);
    vector<int> newScores(atoms.size(), 0);
    vector<int> orderedScores;
    bool goOn = false;
    int n = 0;
    int idx1, idx2;
    int oldTies = atoms.size();
    int newTies = oldTies;
    unsigned int i = 0, j = 0;
    do {
        n++;
        for (i = 0; i < bonds.size(); i++) {
            idx1 = bonds[i]->startAtom->_generalUseN;
            idx2 = bonds[i]->endAtom->_generalUseN;
            newScores[idx1] += oldScores[idx2];
            newScores[idx2] += oldScores[idx1];
        }
        orderedScores = newScores;
        stable_sort(orderedScores.begin(), orderedScores.end());
        newTies = 0;
        for (j = 1; j < orderedScores.size(); j++) {
            if (orderedScores[j] == orderedScores[j - 1])
                newTies++;
        }
        if (newTies < oldTies) {
            goOn = true;
            oldTies = newTies;
            oldScores = newScores;
        } else
            goOn = false;
    } while (goOn);
    return n;
}

CoordgenTemplates sketcherMinimizer::m_templates;
