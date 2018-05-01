/*
 *  sketcherMarchingSquares.h
 *
 *
 *  Created by Nicola Zonta on 13/3/2011.
 *  Copyright Schrodinger, LLC. All rights reserved
 *
 */

#include "sketcherMinimizerMarchingSquares.h"
#include <assert.h>
#include <iostream>
#include "sketcherMinimizer.h"


using namespace std;

sketcherMinimizerMarchingSquares::sketcherMinimizerMarchingSquares()
{
}

sketcherMinimizerMarchingSquares::~sketcherMinimizerMarchingSquares()
{
    clear();
}

void sketcherMinimizerMarchingSquares::initialize(float minx, float maxx,
                                                  float miny, float maxy,
                                                  float x_interval,
                                                  float y_interval)
{
    if (y_interval == 0.f)
        y_interval = x_interval;
    m_xinterval = x_interval;
    m_yinterval = y_interval;

    m_left = minx;
    m_bottom = miny;
    float dx = maxx - minx;
    float dy = maxy - miny;

    assert(dx > 0);
    assert(dy > 0);

    m_XN = (dx / x_interval) + 2;
    m_YN = (dy / y_interval) + 2;
    m_grid.clear();
    m_grid.resize(m_XN * m_YN, 0.f);
    m_lastRowPoints.resize(m_XN, NULL);
}

float sketcherMinimizerMarchingSquares::getNodeValue(unsigned int x,
                                                     unsigned int y) const
{
    if (x + y * m_XN < m_grid.size()) {
        return m_grid[x + y * m_XN];
    } else {
        cerr << "violating grid limits" << endl;
    }
    return 0.f;
}

void sketcherMinimizerMarchingSquares::setValue(float v, unsigned int x,
                                                unsigned int y)
{
    if (x + y * m_XN < m_grid.size()) {
        m_grid[x + y * m_XN] = v;
    } else {
        cerr << "violating grid limits" << endl;
    }
}

void sketcherMinimizerMarchingSquares::clear()
{

    for (unsigned int i = 0; i < m_points.size(); i++)
        delete m_points[i];
    m_points.clear();
    for (unsigned int i = 0; i < m_sides.size(); i++)
        delete m_sides[i];
    m_sides.clear();

    m_grid.clear();
    m_lastRowPoints.clear();
}

void sketcherMinimizerMarchingSquares::setThreshold(float t)
{
    m_threshold = t;
}

float sketcherMinimizerMarchingSquares::getThreshold() const
{
    return m_threshold;
}

float sketcherMinimizerMarchingSquares::toRealx(float x) const
{

    return m_left + (x * m_xinterval);
}

float sketcherMinimizerMarchingSquares::toRealy(float y) const
{

    return m_bottom + (y * m_yinterval);
}

float sketcherMinimizerMarchingSquares::interpolate(float v1, float v2) const
{
    float diff = v2 - v1;
    if (diff < SKETCHER_EPSILON && diff > -SKETCHER_EPSILON)
        return 0.5;
    return ((m_threshold - v1) / (v2 - v1));
}

void sketcherMinimizerMarchingSquares::addSide(
    sketcherMinimizerMarchingSquaresPoint* p1,
    sketcherMinimizerMarchingSquaresPoint* p2)
{
    sketcherMinimizerMarchingSquaresSide* side =
        new sketcherMinimizerMarchingSquaresSide;
    side->p1 = p1;
    side->p2 = p2;
    if (p1->side1)
        p1->side2 = side;
    else
        p1->side1 = side;
    if (p2->side1)
        p2->side2 = side;
    else
        p2->side1 = side;
    m_sides.push_back(side);
}

std::vector<float>
sketcherMinimizerMarchingSquares::getCoordinatesPoints() const
{
    std::vector<float> out;
    for (unsigned int i = 0; i < m_points.size(); i++) {
        out.push_back(m_points[i]->x);
        out.push_back(m_points[i]->y);
    }
    return out;
}

std::vector<std::vector<float>>
sketcherMinimizerMarchingSquares::getOrderedCoordinatesPoints() const
{

    std::vector<std::vector<float>> out;

    bool newShape = true;
    sketcherMinimizerMarchingSquaresPoint* nextPoint = NULL;
    while (newShape) {
        newShape = false;
        nextPoint = NULL;
        for (unsigned int i = 0; i < m_points.size(); i++) {
            if (m_points[i]->visited == false) {
                newShape = true;
                nextPoint = m_points[i];
                break;
            }
        }
        if (nextPoint) {
            std::vector<float> newVec;
            while (nextPoint) {
                nextPoint->visited = true;
                newVec.push_back(nextPoint->x);
                newVec.push_back(nextPoint->y);
                sketcherMinimizerMarchingSquaresPoint* followingPoint1 = NULL;
                if (nextPoint->side1)
                    followingPoint1 = nextPoint->side1->p1;
                if (followingPoint1 == nextPoint)
                    followingPoint1 = nextPoint->side1->p2;

                sketcherMinimizerMarchingSquaresPoint* followingPoint2 = NULL;
                if (nextPoint->side2)
                    followingPoint2 = nextPoint->side2->p1;
                if (followingPoint2 == nextPoint)
                    followingPoint2 = nextPoint->side2->p2;

                bool found = false;
                if (followingPoint1)
                    if (followingPoint1->visited == false) {
                        nextPoint = followingPoint1;
                        found = true;
                    }
                if (!found)
                    if (followingPoint2)
                        if (followingPoint2->visited == false) {
                            nextPoint = followingPoint2;
                            found = true;
                        }

                if (!found)
                    nextPoint = NULL;
            }
            out.push_back(newVec);
        }
    }

    return out;
}

void sketcherMinimizerMarchingSquares::run()
{

    for (unsigned int j = 0; j < m_YN - 1; j++) {
        m_lastCellRightPoint = NULL;

        for (unsigned int i = 0; i < m_XN - 1; i++) {

            assert((i + 1 + j * m_XN) < m_grid.size());
            assert((i + (j + 1) * m_XN) < m_grid.size());
            assert((i + 1 + (j + 1) * m_XN) < m_grid.size());

            // float BL = m_grid [i +     j*m_XN];
            float BR = m_grid[i + 1 + j * m_XN];
            float TL = m_grid[i + (j + 1) * m_XN];
            float TR = m_grid[i + 1 + (j + 1) * m_XN];
            assert(i < m_lastRowPoints.size());
            sketcherMinimizerMarchingSquaresPoint* rp = NULL;
            sketcherMinimizerMarchingSquaresPoint* tp = NULL;
            sketcherMinimizerMarchingSquaresPoint* lp = m_lastCellRightPoint;
            sketcherMinimizerMarchingSquaresPoint* bp = m_lastRowPoints[i];

            if (((BR - m_threshold) * (TR - m_threshold)) < 0) {
                float inter = interpolate(BR, TR);
                float newY = toRealy(inter + j);
                rp = new sketcherMinimizerMarchingSquaresPoint(toRealx(i + 1),
                                                               newY);
                m_points.push_back(rp);
            }
            if (((TL - m_threshold) * (TR - m_threshold)) < 0) {
                float inter = interpolate(TL, TR);
                float newX = toRealx(inter + i);
                tp = new sketcherMinimizerMarchingSquaresPoint(newX,
                                                               toRealy(j + 1));
                m_points.push_back(tp);
            }
            if (rp && tp && lp && bp) {
                if (TL > m_threshold) {
                    addSide(tp, rp);
                    addSide(bp, lp);
                } else {
                    addSide(tp, lp);
                    addSide(bp, rp);
                }
            } else {
                sketcherMinimizerMarchingSquaresPoint *p1 = NULL, *p2 = NULL;
                if (tp)
                    p1 = tp;
                if (rp) {
                    if (p1)
                        p2 = rp;
                    else
                        p1 = rp;
                }
                if (bp) {
                    if (p1)
                        p2 = bp;
                    else
                        p1 = bp;
                }
                if (lp) {
                    if (p1)
                        p2 = lp;
                    else
                        p1 = lp;
                }
                if (p1 && p2)
                    addSide(p1, p2);
            }
            m_lastCellRightPoint = rp;
            m_lastRowPoints[i] = tp;
        }
        m_lastCellRightPoint = NULL;
    }
}