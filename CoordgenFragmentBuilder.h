/*
 Contributors: Nicola Zonta
 Copyright Schrodinger, LLC. All rights reserved
 */

#ifndef COORDGEN_FRAGMENT_BUILDER_H
#define COORDGEN_FRAGMENT_BUILDER_H

#include <iostream>
#include <vector>
#include <stack>
#include <queue>
#include <map>
#include <set>

#include "CoordgenMacrocycleBuilder.h"

class sketcherMinimizerAtom;
class sketcherMinimizerBond;
class sketcherMinimizerRing;
class sketcherMinimizerFragment;
class sketcherMinimizerPointF;

class CoordgenFragmentBuilder
{
  public:
    void initializeCoordinates(sketcherMinimizerFragment* fragment) const;

    bool m_evenAngles; // all bonds are placed at even intervals around the atom
    static  std::vector<sketcherMinimizerAtom*>
    orderRingAtoms(const sketcherMinimizerRing* r);
    static std::vector<sketcherMinimizerAtom*>
    orderChainOfAtoms(const std::vector<sketcherMinimizerAtom*> atoms,
                      sketcherMinimizerAtom* startAtom);
    static std::vector<sketcherMinimizerPointF>
    listOfCoordinatesFromListofRingAtoms(
        const std::vector<sketcherMinimizerAtom*> atoms);
    void setForceOpenMacrocycles(bool b)
    {
        m_macrocycleBuilder.m_forceOpenMacrocycles = b;
    }

    void setPrecision(float f) { m_macrocycleBuilder.setPrecision(f); }

  private:
    std::vector<sketcherMinimizerAtom*>
    getNonSharedAtoms(const std::vector<sketcherMinimizerAtom*> allAtoms,
                      const sketcherMinimizerAtom* pivotAtom1,
                      const sketcherMinimizerAtom* pivotAtom2,
                      const sketcherMinimizerRing* otherRing) const;

    sketcherMinimizerRing* getSharedAtomsWithAlreadyDrawnRing(
        const sketcherMinimizerRing* ring,
        std::vector<sketcherMinimizerAtom*>& fusionAtoms,
        sketcherMinimizerBond*& fusionBond) const;

    void buildRing(sketcherMinimizerRing* ring) const;
    void generateCoordinatesCentralRings(
        std::vector<sketcherMinimizerRing*> centralRings) const;
    sketcherMinimizerRing* findCentralRingOfSystem(
        const std::vector<sketcherMinimizerRing*> rings) const;

    bool findTemplate(const std::vector<sketcherMinimizerRing*> rings) const;
    void generateCoordinatesSideRings(
        std::stack<sketcherMinimizerRing*> sideRings) const;

    void rotateMainFragment(sketcherMinimizerFragment* fragment) const;

    void buildFragment(sketcherMinimizerFragment* fragment) const;
    void buildRings(sketcherMinimizerFragment* fragment) const;
    void buildNonRingAtoms(sketcherMinimizerFragment* fragment) const;
    void
    initializeFusedRingInformation(sketcherMinimizerFragment* fragment) const;
    void
    simplifyRingSystem(const std::vector<sketcherMinimizerRing*> allRings,
                       std::stack<sketcherMinimizerRing*>& sideRings,
                       std::vector<sketcherMinimizerRing*>& centralRings) const;
    void fallbackIfNanCoordinates(sketcherMinimizerFragment* fragment) const;
    void generateCoordinatesNeighborsOfFirstAtomInQueue(
        std::queue<sketcherMinimizerAtom*>& atomQueue,
        std::set<sketcherMinimizerAtom*>& isAtomVisited,
        const sketcherMinimizerFragment* fragment) const;
    std::vector<float>
    neighborsAnglesAtCenter(const sketcherMinimizerAtom* atom) const;
    void initializeVariablesForNeighboursCoordinates(
        sketcherMinimizerAtom* atom,
        std::set<sketcherMinimizerAtom*>& isAtomVisited,
        sketcherMinimizerPointF& startCoordinates,
        std::vector<sketcherMinimizerAtom*>& orderedNeighbours,
        std::vector<float>& angles) const;

    void initializeVariablesForNeighboursCoordinatesRingAtom(
        const sketcherMinimizerAtom* atom,
        std::set<sketcherMinimizerAtom*>& isAtomVisited,
        sketcherMinimizerPointF& startCoordinates,
        std::vector<sketcherMinimizerAtom*>& orderedNeighbours,
        std::vector<float>& angles) const;
    void maybeAddMacrocycleDOF(sketcherMinimizerAtom* atom) const;
    void
    avoidZEInversions(const sketcherMinimizerAtom* at,
                      std::set<sketcherMinimizerAtom*>& isAtomVisited) const;

    float newScorePlanarity(const std::vector<sketcherMinimizerRing*> rings)
        const; // assign a score to the possibility of rings to be drawn on a
               // plane

    CoordgenMacrocycleBuilder m_macrocycleBuilder;
};

#endif /* defined(COORDGEN_FRAGMENT_BUILDER_H) */
