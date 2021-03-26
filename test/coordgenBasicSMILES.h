/*
 * A very, very basic SMILES-like parser. No aromaticity,
 * or zero-order bonds. No chirality/stereochemistry
 *
 * Do not use this as a basis for real SMILES parsers. It just lets
 * us sidestep using a full chemistry toolkit when writing tests.
 */


#include <algorithm>
#include <string>
#include <stack>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>

#include "../sketcherMinimizerMolecule.h"


namespace schrodinger
{

sketcherMinimizerMolecule* approxSmilesParse(const std::string& smiles)
{
    const std::unordered_map<char, int> elements {
        {'H', 1},
        {'C', 6},
        {'N', 7},
        {'S', 16},
        {'O', 8}
    };

    auto mol = new sketcherMinimizerMolecule();
    std::stack<std::stack<sketcherMinimizerAtom*>> tree;
    tree.push({});
    auto* prev = &tree.top();
    std::unordered_map<char, sketcherMinimizerAtom*> cycles;
    int bond_order = 1;

    size_t idx = 0;
    for (auto c: smiles) {
        auto atomic_number = elements.find(c);
        if (atomic_number != elements.end()) {
            auto atom = mol->addNewAtom();
            atom->setAtomicNumber(atomic_number->second);

            if (!prev->empty()) {
                auto bond = mol->addNewBond(atom, prev->top());
                bond->setBondOrder(bond_order);
                bond_order = 1;

            }
            prev->push(atom);
        } else if (c == '=') {
            bond_order = 2;
        } else if (c == '1' || c == '2' || c == '3') {
            auto other = cycles.find(c);
            if (other == cycles.end()) {
                cycles[c] = prev->top();
            } else {
                auto bond = mol->addNewBond(prev->top(), other->second);
                bond->setBondOrder(bond_order);
                bond_order = 1;
                cycles.erase(other);
            }
        } else if (c == '(') {
            auto old = prev->top();
            tree.push({});
            prev = &tree.top();
            prev->push(old);
        } else if (c == ')') {
            tree.pop();
            prev = &tree.top();
        } else {
            std::string msg = "unrecognized symbol: ";
            msg += c;
            msg += " in SMILES: " + smiles;
            throw std::runtime_error(msg);
        }
        ++idx;
    }

    sketcherMinimizerMolecule::assignBondsAndNeighbors(mol->getAtoms(), mol->getBonds());
    return mol;
}

} // namespace schrodinger