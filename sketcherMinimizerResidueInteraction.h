/*
 *  sketcherMinimizerResidueInteraction.h
 *
 *  Created by Nicola Zonta on 13/05/2011.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#ifndef sketcherMINIMIZERRESIDUEINTERACTION_H
#define sketcherMINIMIZERRESIDUEINTERACTION_H

#include "sketcherMinimizerBond.h"

/*class to represent an interaction with protein residues (e.g. h-bond or 
 pi-pi stacking*/
class  sketcherMinimizerResidueInteraction
    : public sketcherMinimizerBond
{
  public:
    sketcherMinimizerResidueInteraction();
    virtual ~sketcherMinimizerResidueInteraction();
    virtual bool isResidueInteraction();

    /*get all the atoms involved at the end side of the interaction*/
    std::vector<sketcherMinimizerAtom*> getAllEndAtoms();
    std::vector<sketcherMinimizerAtom*> getAllStartAtoms();

    /*get all the atoms involved at the start side of the interaction*/
    std::vector<sketcherMinimizerAtom*> m_otherEndAtoms;
    std::vector<sketcherMinimizerAtom*> m_otherStartAtoms;
};

#endif // sketcherMINIMIZERRESIDUEINTERACTION_H
