/*
 *  sketcherMinimizerResidue.cpp
 *
 *  Created by Nicola Zonta on 13/05/2011.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#include "sketcherMinimizerResidue.h"

using namespace std;

sketcherMinimizerResidue::sketcherMinimizerResidue() : sketcherMinimizerAtom()
{
    m_closestLigandAtom = NULL;
}

sketcherMinimizerResidue::~sketcherMinimizerResidue()
{
}

bool sketcherMinimizerResidue::isResidue() const
{
    return true;
}
