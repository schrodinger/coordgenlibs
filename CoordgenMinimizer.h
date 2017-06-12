/*
 Contributors: Nicola Zonta
 Copyright Schrodinger, LLC. All rights reserved
 */

#ifndef COORDGEN_MINIMIZER_H
#define COORDGEN_MINIMIZER_H

#include <iostream>
#include <vector>
#include <set>
#include <map>

class sketcherMinimizerInteraction;
class sketcherMinimizerStretchInteraction;
class sketcherMinimizerBendInteraction;
class sketcherMinimizerClashInteraction;

class sketcherMinimizerMolecule;
class sketcherMinimizerResidue;
class sketcherMinimizerResidueInteraction;
class sketcherMinimizerAtom;
class sketcherMinimizerBond;
class sketcherMinimizerRing;
class sketcherMinimizerFragment;

class CoordgenFragmentDOF;
class CoordgenMinimizer;

class CoordgenDOFSolutions
{
  public:
    CoordgenDOFSolutions(const CoordgenMinimizer* minimizer,
                         sketcherMinimizerMolecule* molecule,
                         std::vector<CoordgenFragmentDOF*> allDofs)
        : m_minimizer(minimizer), m_molecule(molecule), m_allDofs(allDofs)
    {
    }
    float scoreCurrentSolution();
    std::vector<short unsigned int> getCurrentSolution();
    void loadSolution(std::vector<short unsigned int> solution);
    std::pair<std::vector<short unsigned int>, float> findBestSolution() const;
    bool hasSolution(std::vector<short unsigned int> solution);
    std::vector<CoordgenFragmentDOF*> getAllDofs() { return m_allDofs; }
  private:
    const CoordgenMinimizer* m_minimizer;
    sketcherMinimizerMolecule* m_molecule;
    std::map<std::vector<short unsigned int>, float> m_solutions;
    std::vector<CoordgenFragmentDOF*> m_allDofs;
};

class  CoordgenMinimizer
{
  public:
    /*clear all the interactions loaded in the minimizer and free memory*/
    void clearInteractions();

    std::vector<sketcherMinimizerAtom*> _atoms;
    std::vector<sketcherMinimizerBond*> _bonds;
    bool m_evenAngles;
    std::vector<sketcherMinimizerResidue*> _residues;
    std::vector<sketcherMinimizerResidueInteraction*> _residueInteractions;
    std::vector<sketcherMinimizerFragment*> _fragments;
    std::vector<sketcherMinimizerMolecule*> _molecules;

    CoordgenMinimizer();
    ~CoordgenMinimizer();
    void run(); // run a minimization
    void minimizeResidues();
    void minimizeMolecule(sketcherMinimizerMolecule* molecule);
    void minimizeProteinOnlyLID(
        std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains);
    void minimizeAll();
    void setupInteractionsOnlyResidues();
    void setupInteractionsProteinOnly(
        std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains);
    void setupInteractions(bool intrafragmentClashes = false);
    void addInteractionsOfMolecule(sketcherMinimizerMolecule* molecule,
                                   bool intrafragmentClashes = false);

    float scoreInteractions();
    bool findIntermolecularClashes(sketcherMinimizerMolecule* mol1,
                                   sketcherMinimizerMolecule* mol2,
                                   float threshold);
    bool findIntermolecularClashes(std::vector<sketcherMinimizerMolecule*> mols,
                                   float threshold);
    void fixRingsShape();
    float getPrecision() const;
    void setPrecision(float f);

    float
    scoreClashes(sketcherMinimizerMolecule* molecule,
                 bool residueInteractions = false,
                 bool scoreProximityRelationsOnOppositeSides = false) const;
    float scoreCrossBonds(sketcherMinimizerMolecule* molecule,
                          bool residueInteractions = false) const;
    float scoreDofs(sketcherMinimizerMolecule* molecule) const;
    float scoreAtomsInsideRings() const;
    float scoreProximityRelationsOnOppositeSides() const;

    bool avoidClashes();
    bool avoidClashesOfMolecule(
        sketcherMinimizerMolecule* molecule,
        std::vector<sketcherMinimizerInteraction*> extraInteractions =
            std::vector<sketcherMinimizerInteraction*>());
    bool flipFragments(sketcherMinimizerMolecule* molecule, float& clashE);
    bool runLocalSearch(sketcherMinimizerMolecule* molecule,
                        std::vector<CoordgenFragmentDOF*> dofs, int levels,
                        float& clashE, CoordgenDOFSolutions& solutions);

    bool growSolutions(
        std::set<std::vector<short unsigned int>>& allScoredSolutions,
        int& currentTier,
        std::map<std::vector<short unsigned int>, float>& growingSolutions,
        CoordgenDOFSolutions& solutions, float& bestScore);
    bool runSearch(int tier, CoordgenDOFSolutions& solutions);

    std::vector<std::vector<CoordgenFragmentDOF*>>
    buildTuplesOfDofs(std::vector<CoordgenFragmentDOF*> dofs,
                      unsigned int order) const;
    bool runExhaustiveSearch(sketcherMinimizerMolecule* molecule,
                             std::vector<CoordgenFragmentDOF*> dofs,
                             float& clashE, CoordgenDOFSolutions& solutions);
    void runExhaustiveSearchLevel(
        sketcherMinimizerMolecule* molecule,
        std::vector<CoordgenFragmentDOF*>::iterator iterator,
        std::vector<CoordgenFragmentDOF*>& dofs, float& bestResult, bool& abort,
        CoordgenDOFSolutions& solutions);
    bool bondsClash(sketcherMinimizerBond* bond,
                    sketcherMinimizerBond* bond2) const;
    void avoidTerminalClashes(sketcherMinimizerMolecule* molecule,
                              float& clashE);
    static void maybeMinimizeRings(std::vector<sketcherMinimizerRing*> rings);
    static void avoidInternalClashes(sketcherMinimizerFragment* fragment);

    void buildFromFragments(bool firstTime = false) const;
    void buildMoleculeFromFragments(sketcherMinimizerMolecule* molecule,
                                    bool firstTime = false) const;
    /*Apply forces and take a step in the minimization. Returns false if
     * converged, true if not.*/
    bool applyForces(float maxd = 3);

    std::set<sketcherMinimizerAtom*>
    getChetoCs(std::vector<sketcherMinimizerAtom*> allAtoms);

    std::set<sketcherMinimizerAtom*>
    getAminoNs(std::vector<sketcherMinimizerAtom*> allAtoms);

    std::set<sketcherMinimizerAtom*>
    getAlphaCs(std::vector<sketcherMinimizerAtom*> allAtoms,
               std::set<sketcherMinimizerAtom*> chetoCs,
               std::set<sketcherMinimizerAtom*> aminoNs);

    static void checkForClashes(sketcherMinimizerAtom* a);

    bool skipMinimization, skipAvoidClashes, skipFlipFragments,
        m_scoreResidueInteractions;
    static bool hasNaNCoordinates(std::vector<sketcherMinimizerAtom*> atoms);
    bool hasNaNCoordinates();
    static bool
    hasValid3DCoordinates(std::vector<sketcherMinimizerAtom*> atoms);
    static void
    fallbackOn3DCoordinates(std::vector<sketcherMinimizerAtom*> atoms);

    void addExtraInteraction(sketcherMinimizerMolecule* molecule,
                             sketcherMinimizerInteraction* interaction);

  private:
    void addClashInteractionsOfMolecule(sketcherMinimizerMolecule* molecule,
                                        bool intrafragmentClashes);
    void addStretchInteractionsOfMolecule(sketcherMinimizerMolecule* molecule);
    void addBendInteractionsOfMolecule(sketcherMinimizerMolecule* molecule);

    void addChiralInversionConstraintsOfMolecule(
        sketcherMinimizerMolecule* molecule);

    void addPeptideBondInversionConstraintsOfMolecule(
        sketcherMinimizerMolecule* molecule);

    void getFourConsecutiveAtomsThatMatchSequence(
        std::vector<std::vector<sketcherMinimizerAtom*>>&
            consecutiveAtomsGroups,
        std::set<sketcherMinimizerAtom*> firstSet,
        std::set<sketcherMinimizerAtom*> secondSet,
        std::set<sketcherMinimizerAtom*> thirdSet,
        std::set<sketcherMinimizerAtom*> fourthSet) const;

    std::vector<sketcherMinimizerInteraction*> _interactions;
    std::vector<sketcherMinimizerStretchInteraction*> _stretchInteractions;
    std::vector<sketcherMinimizerBendInteraction*> _bendInteractions;

    std::vector<sketcherMinimizerClashInteraction*>
        _intramolecularClashInteractions;
    std::vector<sketcherMinimizerInteraction*> _extraInteractions;
    std::map<sketcherMinimizerMolecule*,
             std::vector<sketcherMinimizerInteraction*>>
        _extraInteractionsOfMolecule;

    float m_maxIterations;
    float m_precision;
};

#endif /* defined(COORDGEN_MINIMIZER_H) */
