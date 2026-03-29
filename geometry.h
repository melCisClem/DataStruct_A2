#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "polygon.h"

constexpr double EPS = 1e-9;

double signed_ring_area(const Ring& ring);
double total_signed_area(const Polygon& polygon);

double cross(const Point& a, const Point& b, const Point& c);
bool almost_equal(double a, double b, double eps = EPS);
bool same_point(const Point& a, const Point& b, double eps = EPS);

bool point_on_segment(const Point& p, const Point& a, const Point& b);
bool segments_intersect(const Point& a, const Point& b,
                        const Point& c, const Point& d);

bool ring_has_minimum_vertices(const Ring& ring, int min_vertices = 3);

#endif