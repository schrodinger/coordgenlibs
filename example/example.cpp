#include <iostream>
#include "../sketcherMinimizer.h"

int main()
{
    sketcherMinimizer minimizer;

    sketcherMinimizerAtom* min_at1 = new sketcherMinimizerAtom();
    sketcherMinimizerAtom* min_at2 = new sketcherMinimizerAtom();
    sketcherMinimizerAtom* min_at3 = new sketcherMinimizerAtom();
    sketcherMinimizerAtom* min_at4 = new sketcherMinimizerAtom();
    sketcherMinimizerAtom* min_at5 = new sketcherMinimizerAtom();
    sketcherMinimizerAtom* min_at6 = new sketcherMinimizerAtom();

    sketcherMinimizerBond* min_bo1 = new sketcherMinimizerBond();
    min_bo1->startAtom = min_at1;
    min_bo1->endAtom = min_at2;

    sketcherMinimizerBond* min_bo2 = new sketcherMinimizerBond();
    min_bo2->startAtom = min_at2;
    min_bo2->endAtom = min_at3;

    sketcherMinimizerBond* min_bo3 = new sketcherMinimizerBond();
    min_bo3->startAtom = min_at3;
    min_bo3->endAtom = min_at4;


    sketcherMinimizerBond* min_bo4 = new sketcherMinimizerBond();
    min_bo4->startAtom = min_at4;
    min_bo4->endAtom = min_at5;

    sketcherMinimizerBond* min_bo5 = new sketcherMinimizerBond();
    min_bo5->startAtom = min_at5;
    min_bo5->endAtom = min_at6;

    sketcherMinimizerBond* min_bo6 = new sketcherMinimizerBond();
    min_bo6->startAtom = min_at6;
    min_bo6->endAtom = min_at1;

    sketcherMinimizerMolecule *min_mol = new sketcherMinimizerMolecule();
    min_mol->_atoms.push_back(min_at1);
    min_mol->_atoms.push_back(min_at2);
    min_mol->_atoms.push_back(min_at3);
    min_mol->_atoms.push_back(min_at4);
    min_mol->_atoms.push_back(min_at5);
    min_mol->_atoms.push_back(min_at6);
    min_mol->_bonds.push_back(min_bo1);
    min_mol->_bonds.push_back(min_bo2);
    min_mol->_bonds.push_back(min_bo3);
    min_mol->_bonds.push_back(min_bo4);
    min_mol->_bonds.push_back(min_bo5);
   // min_mol->_bonds.push_back(min_bo6);


    minimizer.initialize(min_mol);
    minimizer.runGenerateCoordinates();
    for (auto atom : minimizer._atoms)
    {
        std::cerr << atom->coordinates<<std::endl;
    }
    return 0;
}
