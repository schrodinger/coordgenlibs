/*
 *  sketcherMarchingSquares.h
 *
 *
 *  Created by Nicola Zonta on 19/11/2010.
 *  Copyright Schrodinger, LLC. All rights reserved
 *
 */

#ifndef sketcherMINIMIZERMARCHINGSQUARES_H
#define sketcherMINIMIZERMARCHINGSQUARES_H

#include <cstddef>
#include <vector>
#include "CoordgenConfig.hpp"


class sketcherMinimizerPointF;

struct sketcherMinimizerMarchingSquaresPoint;

struct sketcherMinimizerMarchingSquaresSide {
    sketcherMinimizerMarchingSquaresPoint* p1;
    sketcherMinimizerMarchingSquaresPoint* p2;
};

struct sketcherMinimizerMarchingSquaresPoint {
  public:
    sketcherMinimizerMarchingSquaresPoint(float ix, float iy)
    {
        x = ix;
        y = iy;
        side1 = NULL;
        side2 = NULL;
        visited = false;
    }
    float x, y;
    sketcherMinimizerMarchingSquaresSide *side1, *side2;
    bool visited;
};

/*implementation of a marching squares algorithm*/
class  sketcherMinimizerMarchingSquares
{
  public:
    EXPORT_COORDGEN sketcherMinimizerMarchingSquares();
    EXPORT_COORDGEN ~sketcherMinimizerMarchingSquares();
    //  inline void clearGrid ();
    void EXPORT_COORDGEN setValue(float v, unsigned int x, unsigned int y);
    void EXPORT_COORDGEN initialize(float minx, float maxx, float miny, float maxy,
                    float x_interval, float y_interval = 0.f);

    void EXPORT_COORDGEN clear();

    void EXPORT_COORDGEN setThreshold(float t);
    float EXPORT_COORDGEN getThreshold() const;

    float EXPORT_COORDGEN toRealx(float x) const;
    float EXPORT_COORDGEN toRealy(float y) const;

    unsigned int getXN() const { return m_XN; };
    unsigned int getYN() const { return m_YN; };

    void EXPORT_COORDGEN run(); // computes the isovalue points and segments

    std::vector<float>

        /*call after run () is executed, returs the coordinates of all the
           isovalue line points [x1, y1, x2, y2 .. xn, yn] in the order they
           were created*/
        EXPORT_COORDGEN getCoordinatesPoints() const;
    std::vector<std::vector<float>>

        /*call after run () is executed. Returns a vector of isovalue closed
           lines [x1, y1, x2, y2 .. xn, yn]. The points are ordered as they
           appear along the line.*/
        EXPORT_COORDGEN getOrderedCoordinatesPoints() const;

    inline std::vector<float> getRawData() const
    {
        return m_grid;
    }; // returns a vector of all the data set with setValue.

    float EXPORT_COORDGEN getNodeValue(unsigned int i, unsigned int j) const;

  private:
    void addSide(sketcherMinimizerMarchingSquaresPoint* p1,
                 sketcherMinimizerMarchingSquaresPoint* p2);
    float m_xinterval, m_yinterval, m_left, m_bottom;
    std::vector<float> m_grid;
    unsigned int m_XN, m_YN;
    float m_threshold;
    std::vector<sketcherMinimizerMarchingSquaresPoint*> m_lastRowPoints;
    sketcherMinimizerMarchingSquaresPoint* m_lastCellRightPoint;
    float interpolate(float v1, float v2) const;
    std::vector<sketcherMinimizerMarchingSquaresPoint*> m_points;
    std::vector<sketcherMinimizerMarchingSquaresSide*> m_sides;
};

#endif // sketcherMINIMIZERMARCHINGSQUARES_H
