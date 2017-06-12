/*
 Contributors: Nicola Zonta
 Copyright Schrodinger, LLC. All rights reserved
 */

#ifndef COORDGEN_MACROCYCLE_BUILDER_H
#define COORDGEN_MACROCYCLE_BUILDER_H

#include <iostream>
#include <vector>
#include <cstdlib>


class sketcherMinimizerAtom;
class sketcherMinimizerRing;
class sketcherMinimizerBond;
class sketcherMinimizerPointF;

struct doubleBondConstraint {
    bool trans;
    int previousAtom, atom1, atom2, followingAtom;
};

struct ringConstraint {
    ringConstraint(int a, sketcherMinimizerRing* r, bool fo)
    {
        ring = r;
        atom = a;
        forceOutside = fo;
    }
    bool forceOutside;
    int atom;
    sketcherMinimizerRing* ring;
};

struct vertexCoords {
    bool operator!=(const vertexCoords& rhs) const
    {
        if (x != rhs.x)
            return true;
        if (y != rhs.y)
            return true;
        if (z != rhs.z)
            return true;
        return false;
    }
    bool operator<(const vertexCoords& rhs) const
    {
        if (x < rhs.x)
            return true;
        if (y < rhs.y)
            return true;
        if (z < rhs.z)
            return true;
        return false;
    }
    friend const vertexCoords operator+(const vertexCoords& v1,
                                        const vertexCoords& v2)
    {
        return vertexCoords(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
    }
    friend const vertexCoords operator-(const vertexCoords& v1,
                                        const vertexCoords& v2)
    {
        return vertexCoords(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
    }

    bool operator==(const vertexCoords& rhs) const
    {
        if (x == rhs.x && y == rhs.y && z == rhs.z)
            return true;
        return false;
    }

    vertexCoords(int ix, int iy, int iz)
    {
        x = ix;
        y = iy;
        z = iz;
    }
    int x, y, z;

  private:
    friend std::ostream& operator<<(std::ostream& os, const vertexCoords& v);
};

struct pathRestraints {
    std::vector<int> heteroAtoms;
    std::vector<std::pair<int, int>> substitutedAtoms;
};

struct pathConstraints {
    std::vector<doubleBondConstraint> doubleBonds;
    std::vector<ringConstraint> ringConstraints;
    std::vector<int> forceOutside;
};

struct hexCoords {
    hexCoords(int ix, int iy)
    {
        x = ix;
        y = iy;
    }
    bool operator==(const hexCoords& rhs) const
    {
        if (x == rhs.x && y == rhs.y)
            return true;
        return false;
    }
    int x, y;
    int distanceFrom(const hexCoords& origin) const
    {
        int dx = abs(x - origin.x);
        int dy = abs(y - origin.y);
        int dz = abs((-x - y) - (-origin.x - origin.y));
        int max = (dx > dy) ? dx : dy;
        return (dz > max) ? dz : max;
    }
    hexCoords rotate30Degrees()
    {
        int z = -x - y;
        return hexCoords(-z, -x);
    }
    vertexCoords toVertexCoords() const { return vertexCoords(x, y, -x - y); }

  private:
    friend std::ostream& operator<<(std::ostream& os, const hexCoords& h);
};

struct Hex { // unit hexagon used to build polyominoes for macrocycle shapes
    Hex(hexCoords coords) : m_coords(coords) {}
    void setCoords(hexCoords coords) { m_coords = coords; }
    int x() const { return m_coords.x; };
    int y() const { return m_coords.y; };
    int z() const { return -x() - y(); }
    hexCoords coords() const { return m_coords; }
    hexCoords m_coords;
    std::vector<hexCoords> neighbors() const;
    static std::vector<hexCoords> neighboringPositions(hexCoords h);
    vertexCoords followingVertex(vertexCoords v) const;
};

class 
    Polyomino // all the functions assume that the polyomino has no holes.
{
  public:
    Polyomino();
    Polyomino(const Polyomino& p);
    ~Polyomino();
    Polyomino& operator=(const Polyomino& rhv);
    bool isTheSameAs(Polyomino& p) const;
    int size() const;
    void clear();
    void
    markOneVertexAsPentagon(); // to get a path with an odd number of vertices
    std::vector<Hex*> vertexNeighbors(vertexCoords v) const;
    std::vector<hexCoords> freeVertexNeighborPositions(vertexCoords v) const;
    std::vector<vertexCoords> getPath() const;
    Hex* getHex(hexCoords coords) const;
    vertexCoords findOuterVertex() const;
    /*find an outer vertex. Used to start path around polyomino.*/

    int hexagonsAtVertex(
        vertexCoords v) const; // number of Hexs present at that vertex
    void buildWithVerticesN(int totVertices); // round-ish shape
    void buildSkewedBoxShape(int x, int y,
                             bool pentagon = false); // squared shape
    void buildRaggedBoxShape(int x, int y, bool pentagon = false);
    void buildRaggedSmallerBoxShape(
        int x, int y,
        bool pentagon = false); // box alternating rows of length x and x-1
    void buildRaggedBiggerBoxShape(
        int x, int y,
        bool pentagon = false); // box alternating rows of length x and x+1
    std::vector<hexCoords> allFreeNeighbors() const;
    int countNeighbors(hexCoords) const;
    void addHex(hexCoords coords);
    void removeHex(hexCoords coords);
    bool isEquivalentWithout(hexCoords c) const;
    /*does removing this hexagon yield another polyomino with the same number of
     vertices?
     true if the hexagon has 3 neighbors all next to each other
     */

    vertexCoords coordinatesOfSubstituent(vertexCoords pos) const;
    /*give the coordinates of an hypotetical substituent bound to an atom in
     * position pos*/

    std::vector<Hex*> m_list; // holds pointers to all hexagons
    std::vector<vertexCoords> pentagonVertices;

  private:
    void setPentagon(vertexCoords c);
    void resizeGrid(int i) const;
    void reassignHexs() const;
    int getIndexInList(hexCoords coords) const;
    std::vector<Hex*> mutable m_grid; // hold pointers to all positions in the
                                      // grid, NULL pointers for empty positions
    mutable int m_gridSize;
};

class  CoordgenMacrocycleBuilder
{
  public:
    CoordgenMacrocycleBuilder() : m_forceOpenMacrocycles(false){};
    ~CoordgenMacrocycleBuilder(){};
    std::vector<sketcherMinimizerPointF>
    newMacrocycle(sketcherMinimizerRing* r,
                  std::vector<sketcherMinimizerAtom*> atoms) const;
    bool openCycleAndGenerateCoords(sketcherMinimizerRing* ring) const;
    sketcherMinimizerBond* findBondToOpen(sketcherMinimizerRing* ring) const;

    // public to be accessed by tests
    std::vector<doubleBondConstraint>
    getDoubleBondConstraints(std::vector<sketcherMinimizerAtom*>& atoms) const;
    bool m_forceOpenMacrocycles;

    float getPrecision() const;
    void setPrecision(float f);

  private:
    float m_precision;
    std::vector<Polyomino> buildSquaredShapes(int totVertices) const;
    std::vector<Polyomino> removeDuplicates(std::vector<Polyomino>& pols) const;

    std::vector<ringConstraint>
    getRingConstraints(std::vector<sketcherMinimizerAtom*>& atoms) const;
    int getNumberOfChildren(sketcherMinimizerAtom* a,
                            sketcherMinimizerAtom* parent) const;
    pathRestraints getPathRestraints(std::vector<sketcherMinimizerAtom*>& atoms)
        const; // try to avoid things like having eteroatoms or substituents
               // pointing inside of the macrocycle
    pathConstraints
    getPathConstraints(std::vector<sketcherMinimizerAtom*>& atoms)
        const; // a shape is not allowed to break a constraint (i.e. invert
               // stereochemistry on a double bond or put a fused ring inside
               // the macrocycle)
    std::vector<Polyomino>
    listOfEquivalent(Polyomino p) const; // build a list of polyominoes with the
                                         // same number of vertices by removing
                                         // hexagons with 3 neighbors
    std::vector<Polyomino> listOfEquivalents(std::vector<Polyomino> l) const;

    bool checkRingConstraints(std::vector<ringConstraint>& ringConstraints,
                              Polyomino& p, std::vector<vertexCoords>& path,
                              std::vector<int>& neighborNs, int& startI) const;
    bool checkDoubleBoundConstraints(
        std::vector<doubleBondConstraint>& dbConstraints,
        std::vector<vertexCoords>& vertices, int& startI) const;
    int scorePathRestraints(pathRestraints& pr, Polyomino& p,
                            std::vector<vertexCoords>& path,
                            std::vector<int>& neighborNs, int& startI) const;
    bool scorePathConstraints(pathConstraints& pc, Polyomino& p,
                              std::vector<vertexCoords>& path,
                              std::vector<int>& neighborNs, int& startI) const;
    int scorePath(Polyomino& p, std::vector<vertexCoords>& path,
                  std::vector<int>& neighborNs, int& startI,
                  pathConstraints& pc, pathRestraints& pr) const;
    int acceptableShapeScore(int numberOfAtoms) const;
    std::vector<int> getVertexNeighborNs(Polyomino& p,
                                         std::vector<vertexCoords>& path) const;
    int getLowestPeriod(std::vector<int>& neighbors)
        const; // return lowest period after which there is a rotation symmetry

    bool matchPolyominoes(std::vector<Polyomino>& pols, pathConstraints& pc,
                          pathRestraints& pr, int& bestP, int& bestScore,
                          int& bestStart, int& checkedMacrocycles) const;
    bool matchPolyomino(Polyomino& p, pathConstraints& pc, pathRestraints& pr,
                        int& bestStart, int& bestScore) const;
    void writePolyominoCoordinates(std::vector<vertexCoords>& path,
                                   std::vector<sketcherMinimizerAtom*> atoms,
                                   int startI) const;
    sketcherMinimizerPointF coordsOfVertex(vertexCoords& v) const;
};

#endif /* defined(COORDGEN_MACROCYCLE_BUILDER_H) */
