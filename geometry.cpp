#include "geometry.h"

#include <algorithm>
#include <cmath>

double signed_ring_area(const Ring& ring) {
    const std::size_t n = ring.vertices.size();
    if (n < 3) {
        return 0.0;
    }

    double area2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const Point& a = ring.vertices[i];
        const Point& b = ring.vertices[(i + 1) % n];
        area2 += a.x * b.y - b.x * a.y;
    }

    return 0.5 * area2;
}

double total_signed_area(const Polygon& polygon) {
    double total = 0.0;
    for (const Ring& ring : polygon.rings) {
        total += signed_ring_area(ring);
    }
    return total;
}

double cross(const Point& a, const Point& b, const Point& c) {
    // (b-a) x (c-a)
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool almost_equal(double a, double b, double eps) {
    return std::fabs(a - b) <= eps;
}

bool same_point(const Point& a, const Point& b, double eps) {
    return almost_equal(a.x, b.x, eps) && almost_equal(a.y, b.y, eps);
}

bool point_on_segment(const Point& p, const Point& a, const Point& b) {
    if (!almost_equal(cross(a, b, p), 0.0)) {
        return false;
    }

    const double min_x = std::min(a.x, b.x) - EPS;
    const double max_x = std::max(a.x, b.x) + EPS;
    const double min_y = std::min(a.y, b.y) - EPS;
    const double max_y = std::max(a.y, b.y) + EPS;

    return p.x >= min_x && p.x <= max_x &&
           p.y >= min_y && p.y <= max_y;
}

bool segments_intersect(const Point& a, const Point& b,
                        const Point& c, const Point& d) {
    const double ab_c = cross(a, b, c);
    const double ab_d = cross(a, b, d);
    const double cd_a = cross(c, d, a);
    const double cd_b = cross(c, d, b);

    // Proper intersection
    const bool ab_straddles =
    (((ab_c > EPS) && (ab_d < -EPS)) || ((ab_c < -EPS) && (ab_d > EPS)));

const bool cd_straddles =
    (((cd_a > EPS) && (cd_b < -EPS)) || ((cd_a < -EPS) && (cd_b > EPS)));

if (ab_straddles && cd_straddles) {
    return true;
}

    // Collinear / touching cases
    if (almost_equal(ab_c, 0.0) && point_on_segment(c, a, b)) return true;
    if (almost_equal(ab_d, 0.0) && point_on_segment(d, a, b)) return true;
    if (almost_equal(cd_a, 0.0) && point_on_segment(a, c, d)) return true;
    if (almost_equal(cd_b, 0.0) && point_on_segment(b, c, d)) return true;

    return false;
}

bool ring_has_minimum_vertices(const Ring& ring, int min_vertices) {
    return static_cast<int>(ring.vertices.size()) >= min_vertices;
}