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


class  sketcherMinimizerResidueInteraction
    : public sketcherMinimizerBond
{
  public:
    sketcherMinimizerResidueInteraction();
    virtual ~sketcherMinimizerResidueInteraction();
    virtual bool isResidueInteraction();
    std::vector<sketcherMinimizerAtom*> getAllEndAtoms();
    std::vector<sketcherMinimizerAtom*> getAllStartAtoms();

    std::vector<sketcherMinimizerAtom*> m_otherEndAtoms;
    std::vector<sketcherMinimizerAtom*> m_otherStartAtoms;
};

#endif // sketcherMINIMIZERRESIDUEINTERACTION_H
