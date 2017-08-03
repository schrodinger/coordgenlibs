/*
 *  sketcherMinimizerResidue.h
 *
 *  Created by Nicola Zonta on 13/05/2011.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#ifndef sketcherMINIMIZERRESIDUE_H
#define sketcherMINIMIZERRESIDUE_H

#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBond.h"


/*class to represent protein residues*/
class  sketcherMinimizerResidue : public sketcherMinimizerAtom
{
  public:
    sketcherMinimizerResidue();
    virtual ~sketcherMinimizerResidue();
    virtual bool isResidue() const;

    /*compute coordinates based on the position of the closest ligand atom*/
    sketcherMinimizerPointF computeStartingCoordinates(float d = 2.f)
    {
        sketcherMinimizerPointF out = templateCoordinates;
        if (m_closestLigandAtom)
            out = m_closestLigandAtom->getSingleAdditionVector() * d +
                  m_closestLigandAtom->coordinates;
        if (residueInteractions.size()) {
            int nn = 0;
            sketcherMinimizerPointF coords(0.f, 0.f);
            for (unsigned int i = 0; i < residueInteractions.size(); i++) {
                sketcherMinimizerAtom* n = residueInteractions[i]->endAtom;
                if (n == this)
                    n = residueInteractions[i]->startAtom;
                if (!n->isResidue()) {
                    coords += n->getSingleAdditionVector() * d + n->coordinates;
                    nn++;
                } else if (n->coordinatesSet) {
                    coords += n->coordinates;
                    nn++;
                }
            }
            if (nn > 0)
                coords /= float(nn);
            out = coords;
        }
        return out;
    }

    std::string chain;
    int resnum;
    sketcherMinimizerAtom* m_closestLigandAtom;
};

#endif // sketcherMINIMIZERRESIDUE_H
