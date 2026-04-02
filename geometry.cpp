#include "geometry.h"

#include <algorithm>
#include <cmath>

/**
 * Calculates the signed area of a polygonal ring using the Shoelace Formula.
 * The area is calculated by summing the cross products of sequential vertex vectors.
 * Counter-clockwise (CCW) orientation returns a positive area.
 * Clockwise (CW) orientation returns a negative area.
 */
double signed_ring_area(const Ring& ring) {
    const std::size_t n = ring.vertices.size();
    if (n < 3) {
        return 0.0;
    }

    double area2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const Point& a = ring.vertices[i];
        const Point& b = ring.vertices[(i + 1) % n];
        // Cross product sum (Shoelace formula)
        area2 += a.x * b.y - b.x * a.y;
    }

    return 0.5 * area2;
}

/**
 * Calculates the total signed area of a polygon by summing the signed areas of all its rings.
 * Following GIS conventions, the exterior ring is usually CCW (+) and holes are CW (-).
 * The sum represents the net "filled" area of the polygon.
 */
double total_signed_area(const Polygon& polygon) {
    double total = 0.0;
    for (const Ring& ring : polygon.rings) {
        total += signed_ring_area(ring);
    }
    return total;
}

/**
 * Computes the 2D cross product of vectors (b-a) and (c-a).
 * The sign of the result indicates the orientation of point c relative to line ab:
 * Positive: c is to the left of ab (CCW).
 * Negative: c is to the right of ab (CW).
 * Zero: a, b, and c are collinear.
 */
double cross(const Point& a, const Point& b, const Point& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

/**
 * Checks if two doubles are equal within a small floating-point epsilon.
 * Handles precision issues inherent in floating-point arithmetic.
 */
bool almost_equal(double a, double b, double eps) {
    return std::fabs(a - b) <= eps;
}

/**
 * Checks if two points are identical within a small epsilon.
 * Two points are considered the same if both their x and y coordinates are almost equal.
 */
bool same_point(const Point& a, const Point& b, double eps) {
    return almost_equal(a.x, b.x, eps) && almost_equal(a.y, b.y, eps);
}

/**
 * Checks if point p lies on the line segment defined by endpoints a and b.
 * The point must be collinear with ab and fall within the bounding box of the segment.
 */
bool point_on_segment(const Point& p, const Point& a, const Point& b) {
    // Check for collinearity using the cross product
    if (!almost_equal(cross(a, b, p), 0.0)) {
        return false;
    }

    // Check if the point lies within the bounding box of the segment
    const double min_x = std::min(a.x, b.x) - EPS;
    const double max_x = std::max(a.x, b.x) + EPS;
    const double min_y = std::min(a.y, b.y) - EPS;
    const double max_y = std::max(a.y, b.y) + EPS;

    return p.x >= min_x && p.x <= max_x &&
           p.y >= min_y && p.y <= max_y;
}

/**
 * Determines if line segments (a,b) and (c,d) intersect.
 * Employs the "straddle test": segments intersect if their endpoints lie on
 * opposite sides of the other segment's line.
 * Also handles special cases where an endpoint lies exactly on the other segment.
 */
bool segments_intersect(const Point& a, const Point& b,
                        const Point& c, const Point& d) {
    const double ab_c = cross(a, b, c);
    const double ab_d = cross(a, b, d);
    const double cd_a = cross(c, d, a);
    const double cd_b = cross(c, d, b);

    // Proper intersection: endpoints strictly straddle the other line
    const bool ab_straddles =
        (((ab_c > EPS) && (ab_d < -EPS)) || ((ab_c < -EPS) && (ab_d > EPS)));

    const bool cd_straddles =
        (((cd_a > EPS) && (cd_b < -EPS)) || ((cd_a < -EPS) && (cd_b > EPS)));

    if (ab_straddles && cd_straddles) {
        return true;
    }

    // Boundary cases: one segment's endpoint lies on the other segment
    if (almost_equal(ab_c, 0.0) && point_on_segment(c, a, b)) return true;
    if (almost_equal(ab_d, 0.0) && point_on_segment(d, a, b)) return true;
    if (almost_equal(cd_a, 0.0) && point_on_segment(a, c, d)) return true;
    if (almost_equal(cd_b, 0.0) && point_on_segment(b, c, d)) return true;

    return false;
}

/**
 * Checks if a ring has at least the minimum number of vertices to remain a valid polygon.
 * A simple polygon must have at least 3 vertices to enclose an area.
 */
bool ring_has_minimum_vertices(const Ring& ring, int min_vertices) {
    return static_cast<int>(ring.vertices.size()) >= min_vertices;
}
