/*
 *  sketcherMinimizer.h
 *
 *  Created by Nicola Zonta on 13/04/2010.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */
#ifndef sketcherMINIMIZER
#define sketcherMINIMIZER

#include <map>
#include <vector>
#include <stack>

#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerResidue.h"
#include "sketcherMinimizerResidueInteraction.h"

#include "sketcherMinimizerStretchInteraction.h"
#include "sketcherMinimizerClashInteraction.h"
#include "sketcherMinimizerBendInteraction.h"

#include "sketcherMinimizerFragment.h"
#include "sketcherMinimizerMolecule.h"
#include "sketcherMinimizerRing.h"
#include "sketcherMinimizerMarchingSquares.h"
#include <iostream>


#include "CoordgenFragmentBuilder.h"
#include "CoordgenMinimizer.h"

static const float SKETCHER_STANDARD_PRECISION = 1.f;
static const float SKETCHER_QUICK_PRECISION = 0.2;
static const float SKETCHER_BEST_PRECISION = 3;

class sketcherAtom;
class sketcherBond;
class sketcherMolecule;
class sketcherMinimimizerInteraction;

typedef struct {
    std::vector<sketcherMinimizerPointF> additionVectors;
    std::vector<sketcherMinimizerPointF> centers;
    std::vector<int> counters;
} proximityData;

class  CoordgenTemplates
{
  public:
    CoordgenTemplates() {}
    ~CoordgenTemplates()
    {
        for (auto molecule : m_templates) {
            for (auto atom : molecule->_atoms) {
                delete atom;
            }
            for (auto bond : molecule->_bonds) {
                delete bond;
            }
            delete molecule;
        }
        m_templates.clear();
    }
    std::vector<sketcherMinimizerMolecule*>& getTemplates()
    {
        return m_templates;
    }

  private:
    std::vector<sketcherMinimizerMolecule*> m_templates;
};

class  sketcherMinimizer
{

  public:
    sketcherMinimizer(float precision = SKETCHER_STANDARD_PRECISION);
    ~sketcherMinimizer();

    CoordgenFragmentBuilder m_fragmentBuilder;
    CoordgenMinimizer m_minimizer;


   // bool generateCoordinates(ChmMol& molecule);

    bool runGenerateCoordinates(); // run coordinates generation and return true
                                   // if the pose is considered optimal

    bool structurePassSanityCheck() const;
    void clear();
    void initialize(sketcherMinimizerMolecule* minMol);
    static void canonicalOrdering(sketcherMinimizerMolecule* minMol);

   // void initializeFromMolecule(ChmMol& mol);
    void splitIntoMolecules(sketcherMinimizerMolecule* mol,
                            std::vector<sketcherMinimizerMolecule*>& mols);
    void flagCrossAtoms();

    void minimizeAll();
    void minimizeMolecule(sketcherMinimizerMolecule* molecule);

    void bestRotation();
    void
    addBestRotationInfoForPeptides(std::vector<std::pair<float, float>>& angles,
                                   std::vector<sketcherMinimizerAtom*> atoms);

    void maybeFlipPeptides(std::vector<sketcherMinimizerAtom*> atoms,
                           float& scoreX);
    void maybeFlip();
    void assignPseudoZ();
    void writeStereoChemistry();
    void arrangeMultipleMolecules();
    void placeMoleculesWithProximityRelations(
        std::vector<sketcherMinimizerMolecule*> proximityMols);
    void placeMolResidueLigandStyle(sketcherMinimizerMolecule* mol,
                                    sketcherMinimizerMolecule* parent);
    void flipIfCrossingInteractions(sketcherMinimizerMolecule* mol);

    std::vector<proximityData> buildProximityDataVector(
        std::vector<sketcherMinimizerMolecule*>& proximityMols,
        std::map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap);
    void translateMoleculesWithProximityRelations(
        std::vector<sketcherMinimizerMolecule*>& proximityMols,
        std::map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap,
        std::map<sketcherMinimizerMolecule*, sketcherMinimizerPointF>&
            templateCenters,
        std::vector<proximityData>& proximityDataVecto);
    void rotateMoleculesWithProximityRelations(
        std::vector<sketcherMinimizerMolecule*>& proximityMols,
        std::map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap,
        std::vector<proximityData>& proximityDataVector);

    void placeResiduesInCrowns();

    bool fillShape(std::vector<std::vector<sketcherMinimizerResidue*>>& SSEs,
                   const std::vector<sketcherMinimizerPointF>& shape,
                   int shapeN);

    void placeSSE(std::vector<sketcherMinimizerResidue*> SSE,
                  const std::vector<sketcherMinimizerPointF>& shape, int shapeN,
                  std::vector<bool>& penalties,
                  std::set<sketcherMinimizerResidue*>& outliers,
                  bool placeOnlyInteracting = false);

    float scoreSSEPosition(std::vector<sketcherMinimizerResidue*> SSE,
                           const std::vector<sketcherMinimizerPointF>& shape,
                           int shapeN, std::vector<bool>& penalties, float f,
                           float increment);

    float scoreSSEBondStretch(sketcherMinimizerPointF coordinates1,
                              sketcherMinimizerPointF coordinates2);

    float getResidueDistance(float startF, float increment,
                             sketcherMinimizerResidue* res,
                             std::vector<sketcherMinimizerResidue*> SSE);
    int getShapeIndex(std::vector<sketcherMinimizerPointF> shape, float f);

    void markSolution(std::pair<float, float> solution,
                      std::vector<sketcherMinimizerResidue*> SSE,
                      const std::vector<sketcherMinimizerPointF>& shape,
                      std::vector<bool>& penalties,
                      std::set<sketcherMinimizerResidue*>& outliers);

    std::vector<sketcherMinimizerPointF> shapeAroundLigand(int crownN);
    std::vector<std::vector<sketcherMinimizerResidue*>>
    groupResiduesInSSEs(std::vector<sketcherMinimizerResidue*> residues);

    float
    scoreResiduePosition(int index,
                         const std::vector<sketcherMinimizerPointF>& shape,
                         int shapeN, std::vector<bool>& penalties,
                         sketcherMinimizerResidue* residue);

    void placeResidues(std::vector<sketcherMinimizerAtom*> atoms =
                           std::vector<sketcherMinimizerAtom*>(0));
    void placeResiduesProteinOnlyMode();

    void placeResiduesProteinOnlyModeCircleStyle(
        std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains);

    void placeResiduesProteinOnlyModeLIDStyle(
        std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains);

    std::vector<sketcherMinimizerResidue*> orderResiduesOfChains(
        std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains);

    std::map<std::string, sketcherMinimizerPointF>
    computeChainsStartingPositionsMetaMol(
        std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains);
    void shortenInteractions(
        std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains);
    sketcherMinimizerPointF exploreGridAround(
        sketcherMinimizerPointF centerOfGrid, unsigned int levels, float gridD,
        float dx = 0.f, float dy = 0.f, float distanceFromAtoms = -1.f,
        bool watermap = false,
        sketcherMinimizerResidue* residueForInteractions = NULL,
        sketcherMinimizerPointF direction = sketcherMinimizerPointF(0, 1));
    sketcherMinimizerPointF exploreMolPosition(sketcherMinimizerMolecule* mol,
                                               unsigned int levels, float gridD,
                                               float distanceFromAtoms = -1.f);

    void addToVector(float weight, float angle,
                     std::vector<std::pair<float, float>>& angles);

    static float
    testAlignment(sketcherMinimizerPointF direction,
                  std::pair<sketcherMinimizerPointF, float> templat);
    static sketcherMinimizerPointF scoreDirections(
        sketcherMinimizerFragment* fragment, float angle,
        std::vector<std::pair<sketcherMinimizerPointF, float>> directions,
        bool& invert);
    static void alignWithParentDirection(sketcherMinimizerFragment* f,
                                         sketcherMinimizerPointF position,
                                         float angle);
    static bool
    alignWithParentDirectionConstrained(sketcherMinimizerFragment* fragment,
                                        sketcherMinimizerPointF position,
                                        float angle);
    static bool
    alignWithParentDirectionUnconstrained(sketcherMinimizerFragment* fragment,
                                          float angle);
    static std::vector<sketcherMinimizerBond*>
    getAllTerminalBonds(sketcherMinimizerFragment* fragment);
    static std::vector<std::pair<sketcherMinimizerPointF, float>>
    findDirectionsToAlignWith(sketcherMinimizerFragment* fragment);

    std::vector<sketcherMinimizerAtom*> _atoms;
    std::vector<sketcherMinimizerAtom*> _referenceAtoms;
    std::vector<sketcherMinimizerResidue*> _residues;
    std::vector<sketcherMinimizerResidueInteraction*> _residueInteractions;

    std::vector<sketcherMinimizerFragment*> _fragments;
    std::vector<sketcherMinimizerFragment*> _independentFragments;

    std::vector<sketcherMinimizerBond*> _bonds;
    std::vector<sketcherMinimizerBond*> _referenceBonds;
    std::vector<sketcherMinimizerBond*> m_proximityRelations;
    std::vector<sketcherMinimizerBond*> m_extraBonds;
    std::vector<sketcherMinimizerMolecule*> _molecules;

    void assignLongestChainFromHere(sketcherMinimizerFragment* f);
    void assignNumberOfChildrenAtomsFromHere(sketcherMinimizerFragment* f);

//    void exportCoordinates(ChmMol& molecule);

    void findFragments();
    void initializeFragments();

    void constrainAllAtoms();
    void constrainAtoms(std::vector<bool> constrained);
    void fixAtoms(std::vector<bool> fixed);

    void setScoreResidueInteractions(bool b);

    static sketcherMinimizerAtom*
    pickBestAtom(std::vector<sketcherMinimizerAtom*>& atoms);

    static sketcherMinimizerRing* sameRing(const sketcherMinimizerAtom* at1,
                                           const sketcherMinimizerAtom* at2,
                                           const sketcherMinimizerAtom* at3);
    static sketcherMinimizerRing* sameRing(const sketcherMinimizerAtom* at1,
                                           const sketcherMinimizerAtom* at2);

    static sketcherMinimizerBond* getBond(const sketcherMinimizerAtom* a1,
                                          const sketcherMinimizerAtom* a2);

    void findClosestAtomToResidues(std::vector<sketcherMinimizerAtom*> atoms =
                                       std::vector<sketcherMinimizerAtom*>(0));

    static float RMSD(std::vector<sketcherMinimizerPointF> templates,
                      std::vector<sketcherMinimizerPointF> points);

    static void svd(float* a, float* U, float* Sig,
                    float* V); // singular value decomposition for 2x2 matrices.
                               // used for 2D alignment.
    static void alignmentMatrix(std::vector<sketcherMinimizerPointF> ref,
                                std::vector<sketcherMinimizerPointF> points,
                                float* m);

    static void
    checkIdentity(std::vector<unsigned int> solution, int newSol,
                  std::vector<bool>& matrix,
                  std::vector<sketcherMinimizerPointF>& templateCoordinates,
                  std::vector<std::vector<int>>& molBonds,
                  std::vector<std::vector<int>>& templateBonds,
                  std::vector<std::vector<int>>& molCisTransChains,
                  std::vector<bool>& molIsCis, unsigned int size, bool& found,
                  std::vector<unsigned int>& mapping);
    static bool compare(std::vector<sketcherMinimizerAtom*> atoms,
                        std::vector<sketcherMinimizerBond*> bonds,
                        sketcherMinimizerMolecule* templ,
                        std::vector<unsigned int>& mapping);
    static int morganScores(std::vector<sketcherMinimizerAtom*> atoms,
                            std::vector<sketcherMinimizerBond*> bonds,
                            std::vector<int>& scores);

    std::string m_chainHint;

    static void loadTemplates();
    static CoordgenTemplates m_templates;
};

#endif // sketcherMINIMIZER
