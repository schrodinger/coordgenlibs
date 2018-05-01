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
#include "CoordgenConfig.hpp"

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


/*class to handle templates of common difficult ring structures*/
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
    void setTemplateDir(std::string dir)
    {
        m_templateDir = dir;
    }
    std::string getTemplateDir()
    {
        return m_templateDir;
    }


  private:
    std::vector<sketcherMinimizerMolecule*> m_templates;
    std::string m_templateDir = "";
};


/*main class. Creates 2d coordinates for molecular inputs*/
class  sketcherMinimizer
{

  public:
    EXPORT_COORDGEN sketcherMinimizer(float precision = SKETCHER_STANDARD_PRECISION);
    EXPORT_COORDGEN ~sketcherMinimizer();

    CoordgenFragmentBuilder m_fragmentBuilder;
    CoordgenMinimizer m_minimizer;

    /*run coordinates generation and return true
      if the pose is considered optimal*/
    bool EXPORT_COORDGEN runGenerateCoordinates();

    /*
     return true if the molecules structure is reasonable (e.g. reasonable amount of fused rings)
     */
    bool structurePassSanityCheck() const;

    /*
     clear data and free memory
     */
    void clear();

    /*
     initialize data from given molecule
     */
    void EXPORT_COORDGEN initialize(sketcherMinimizerMolecule* minMol);

    /*put atoms in a canonical order to reduce dependency from order in the input vector */
    static void canonicalOrdering(sketcherMinimizerMolecule* minMol);

   // void initializeFromMolecule(ChmMol& mol);

    /*if mol contains separate molecules, split them into a vector*/
    void splitIntoMolecules(sketcherMinimizerMolecule* mol,
                            std::vector<sketcherMinimizerMolecule*>& mols);

    /*flag atoms that will be drawn with 90° angles (e.g. phosphate P)*/
    void flagCrossAtoms();

    /*assign coordinates to all molecules and residues*/
    void minimizeAll();


    /*assign coordinates to given molecule*/
    void minimizeMolecule(sketcherMinimizerMolecule* molecule);

    /*find the best angle to rotate each molecule*/
    void bestRotation();

    /*add info to choose the best angle so that, if present, peptide chains are horizontal*/
    void
    addBestRotationInfoForPeptides(std::vector<std::pair<float, float>>& angles,
                                   std::vector<sketcherMinimizerAtom*> atoms);


    /*if a peptide chain is present make sure that the N term is on the left*/
    void maybeFlipPeptides(std::vector<sketcherMinimizerAtom*> atoms,
                           float& scoreX);

    /*consider flipping the molecule horizontaly and/or vertically*/
    void maybeFlip();

    /*mark atoms with wedges as above or below the plane to correctly draw crossing bonds*/
    void assignPseudoZ();

    /*write wedges and dashed bonds to mark stereochemistry*/
    void EXPORT_COORDGEN writeStereoChemistry();

    /*arrange multiple molecules next to each other*/
    void arrangeMultipleMolecules();

    /*arrange molecules that have parts that interact with each other so that they are close*/
    void placeMoleculesWithProximityRelations(
        std::vector<sketcherMinimizerMolecule*> proximityMols);

    /*place molecules so that one is in the middle and other are around*/
    void placeMolResidueLigandStyle(sketcherMinimizerMolecule* mol,
                                    sketcherMinimizerMolecule* parent);

    /* if the molecule has more than one interaction and they cross mirror its
     coordinates so they don't cross anymore*/
    void flipIfCrossingInteractions(sketcherMinimizerMolecule* mol);

    /*build data vectors to place molecules with proximity relations*/
    std::vector<proximityData> buildProximityDataVector(
        std::vector<sketcherMinimizerMolecule*>& proximityMols,
        std::map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap);

    /*translate molecules with proximity relations*/
    void translateMoleculesWithProximityRelations(
        std::vector<sketcherMinimizerMolecule*>& proximityMols,
        std::map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap,
        std::map<sketcherMinimizerMolecule*, sketcherMinimizerPointF>&
            templateCenters,
        std::vector<proximityData>& proximityDataVecto);

    /*rotate molecules with proximity relations*/
    void rotateMoleculesWithProximityRelations(
        std::vector<sketcherMinimizerMolecule*>& proximityMols,
        std::map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap,
        std::vector<proximityData>& proximityDataVector);

    /*place residues in concentric contours around the ligand*/
    void placeResiduesInCrowns();

    /*place residues from the given strands of consecutive residues to fill
     the given path*/
    bool fillShape(std::vector<std::vector<sketcherMinimizerResidue*>>& SSEs,
                   const std::vector<sketcherMinimizerPointF>& shape,
                   int shapeN);

    /*place a single strand of consecutive resiudes*/
    void placeSSE(std::vector<sketcherMinimizerResidue*> SSE,
                  const std::vector<sketcherMinimizerPointF>& shape, int shapeN,
                  std::vector<bool>& penalties,
                  std::set<sketcherMinimizerResidue*>& outliers,
                  bool placeOnlyInteracting = false);

    /*score the position of the given strands*/
    float scoreSSEPosition(std::vector<sketcherMinimizerResidue*> SSE,
                           const std::vector<sketcherMinimizerPointF>& shape,
                           int shapeN, std::vector<bool>& penalties, float f,
                           float increment);

    /*score the distance between the two given points of connected residues*/
    float scoreSSEBondStretch(sketcherMinimizerPointF coordinates1,
                              sketcherMinimizerPointF coordinates2);

    /*return the position of res, which is part of SSE, given that the
     first residue of SSE is placed at startF and consecutive residues are placed
     increment away from each other. All distances are expressed in floats, where
     0.f is
     an arbitrary starting point, 0.5 is the opposite side of the curve and 1.0 is
     again
     the starting point*/
    float getResidueDistance(float startF, float increment,
                             sketcherMinimizerResidue* res,
                             std::vector<sketcherMinimizerResidue*> SSE);

    /*return the vector index corresponding to floatPosition*/
    int getShapeIndex(std::vector<sketcherMinimizerPointF> shape, float floatPosition);

    /*solution represent the placement chosen for residues in SSE. Mark the
     corresponding
     sections of the crown to prevent other residues to be placed there*/
    void markSolution(std::pair<float, float> solution,
                      std::vector<sketcherMinimizerResidue*> SSE,
                      const std::vector<sketcherMinimizerPointF>& shape,
                      std::vector<bool>& penalties,
                      std::set<sketcherMinimizerResidue*>& outliers);

    /*return a concentric shape around the ligand. CrownN controls how far away
     from the ligand the shape is*/
    std::vector<sketcherMinimizerPointF> shapeAroundLigand(int crownN);

    /*group residues in strands of consecutive residues*/
    std::vector<std::vector<sketcherMinimizerResidue*>>
    groupResiduesInSSEs(std::vector<sketcherMinimizerResidue*> residues);

    /*score the position of given residues*/
    float
    scoreResiduePosition(int index,
                         const std::vector<sketcherMinimizerPointF>& shape,
                         int shapeN, std::vector<bool>& penalties,
                         sketcherMinimizerResidue* residue);

    /*assign coordinates to residues*/
    void placeResidues(std::vector<sketcherMinimizerAtom*> atoms =
                           std::vector<sketcherMinimizerAtom*>(0));

    /*assign coordinates to residues in the context of a protein-protein
     interaction diagram*/
    void placeResiduesProteinOnlyMode();


    /*assign coordinates to residues in a protein-protein interaction
     diagram shaped as a circle*/
    void placeResiduesProteinOnlyModeCircleStyle(
        std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains);


    /*assign coordinates to residues in a protein-protein interaction
     diagram shaped as a LID*/
    void placeResiduesProteinOnlyModeLIDStyle(
        std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains);

    /*
     order residues for drawing, so that residues interacting together are drawn one
     after the other and
     residues with more interactions are drawn first
     */
    std::vector<sketcherMinimizerResidue*> orderResiduesOfChains(
        std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains);


    /*find center position for each chain of residues using a meta-molecule approach,
     building a molecule where each atom represents a chain and each bond connects
     two interacting chains*/
    std::map<std::string, sketcherMinimizerPointF>
    computeChainsStartingPositionsMetaMol(
        std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains);

    /*place interacting residues closer to each other, so they end up at the perifery
     of the chain*/
    void shortenInteractions(
        std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains);

    /*explore positions in a grid around the current one to ease clashes*/
    sketcherMinimizerPointF exploreGridAround(
        sketcherMinimizerPointF centerOfGrid, unsigned int levels, float gridD,
        float dx = 0.f, float dy = 0.f, float distanceFromAtoms = -1.f,
        bool watermap = false,
        sketcherMinimizerResidue* residueForInteractions = NULL,
        sketcherMinimizerPointF direction = sketcherMinimizerPointF(0, 1));

    sketcherMinimizerPointF exploreMolPosition(sketcherMinimizerMolecule* mol,
                                               unsigned int levels, float gridD,
                                               float distanceFromAtoms = -1.f);


    /*add given angle to a vector of angles, clustering them
     if a close angle is already in the vector*/
    void addToVector(float weight, float angle,
                     std::vector<std::pair<float, float>>& angles);


    /*return a score of the alignment between direction and templat.first, weight on the
     angle between the two and templat.second*/
    static float
    testAlignment(sketcherMinimizerPointF direction,
                  std::pair<sketcherMinimizerPointF, float> templat);


    /*find the best alignment of a fragment to its parent and set invert in case the fragment
     needs to be flipped*/
    static sketcherMinimizerPointF scoreDirections(
        sketcherMinimizerFragment* fragment, float angle,
        std::vector<std::pair<sketcherMinimizerPointF, float>> directions,
        bool& invert);

    /*align the fragment to its parent*/
    static void alignWithParentDirection(sketcherMinimizerFragment* f,
                                         sketcherMinimizerPointF position,
                                         float angle);

    /*align the fragment to its parent in the case of constrained coordinates*/
    static bool
    alignWithParentDirectionConstrained(sketcherMinimizerFragment* fragment,
                                        sketcherMinimizerPointF position,
                                        float angle);

    /*align the fragment to its parent in the case of unconstrained coordinates*/
    static bool
    alignWithParentDirectionUnconstrained(sketcherMinimizerFragment* fragment,
                                          float angle);

    /*get all bonds to a terminal atom*/
    static std::vector<sketcherMinimizerBond*>
    getAllTerminalBonds(sketcherMinimizerFragment* fragment);

    /*return a list of vectors the given fragment can be aligned with and a score of 
     the importance of each*/
    static std::vector<std::pair<sketcherMinimizerPointF, float>>
    findDirectionsToAlignWith(sketcherMinimizerFragment* fragment);

    std::vector<sketcherMinimizerAtom*> getAtoms() {return _atoms;}

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

    /*split molecules into rigid fragments*/
    void EXPORT_COORDGEN findFragments();

    /*initialize data and coordinates for each fragment*/
    void initializeFragments();

    /*constrain coordinates on all atoms*/
    void EXPORT_COORDGEN constrainAllAtoms();

    /*constrain coordinates on atoms corresponding to true*/
    void EXPORT_COORDGEN constrainAtoms(std::vector<bool> constrained);

    /*fix cooordinates (i.e. guarantee they will not change) on atoms marked as true*/
    void EXPORT_COORDGEN fixAtoms(std::vector<bool> fixed);

    /*set a flag to enable/disable the scoring of interactions with residues*/
    void EXPORT_COORDGEN setScoreResidueInteractions(bool b);

    /*
     pick one atom out of the vector. Arbitrary criteria such as atomic number and connectivity
     are used
     */
    static sketcherMinimizerAtom*
    pickBestAtom(std::vector<sketcherMinimizerAtom*>& atoms);


    /*if the three atoms share a ring, return it*/
    static sketcherMinimizerRing* sameRing(const sketcherMinimizerAtom* at1,
                                           const sketcherMinimizerAtom* at2,
                                           const sketcherMinimizerAtom* at3);

    /*if the two atoms share a ring, return it*/
    static sketcherMinimizerRing* sameRing(const sketcherMinimizerAtom* at1,
                                           const sketcherMinimizerAtom* at2);

    /*if the two atoms share a bond, return it*/
    static sketcherMinimizerBond* getBond(const sketcherMinimizerAtom* a1,
                                          const sketcherMinimizerAtom* a2);

    /*for each residue, find the closest atom among atoms, or all atoms
     if none are given*/
    void findClosestAtomToResidues(std::vector<sketcherMinimizerAtom*> atoms =
                                       std::vector<sketcherMinimizerAtom*>(0));

    /*calculate root mean square deviation between templates and points*/
    static float RMSD(std::vector<sketcherMinimizerPointF> templates,
                      std::vector<sketcherMinimizerPointF> points);

    /* singular value decomposition for 2x2 matrices.
     used for 2D alignment.*/
    static void svd(float* a, float* U, float* Sig,
                    float* V);


    /*set m to a rotation matrix to align ref to points*/
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


    /*compare atoms and bonds to template and map which atom is which in case of a positive match*/
    static bool compare(std::vector<sketcherMinimizerAtom*> atoms,
                        std::vector<sketcherMinimizerBond*> bonds,
                        sketcherMinimizerMolecule* templ,
                        std::vector<unsigned int>& mapping);

    /*calculate morgan scores for the given input*/
    static int morganScores(std::vector<sketcherMinimizerAtom*> atoms,
                            std::vector<sketcherMinimizerBond*> bonds,
                            std::vector<int>& scores);

    std::string m_chainHint;

    /*load the templates from the template file*/
    static void EXPORT_COORDGEN setTemplateFileDir(std::string dir);
    static void EXPORT_COORDGEN loadTemplates();
    static EXPORT_COORDGEN CoordgenTemplates m_templates;
};

#endif // sketcherMINIMIZER
