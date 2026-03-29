#ifndef SIMPLIFIER_H
#define SIMPLIFIER_H

#include "polygon.h"

struct CollapseCandidate {
    bool valid{false};
    Point e{};
    double areal_displacement{0.0};
    int start_index{-1};
};

Polygon simplify_polygon(const Polygon& input, int target_vertices, double& areal_displacement);

double triangle_signed_area(const Point& a, const Point& b, const Point& c);
double triangle_area_abs(const Point& a, const Point& b, const Point& c);
CollapseCandidate compute_collapse_candidate(const Ring& ring, int start_index);
void apply_collapse(Ring& ring, const CollapseCandidate& c);
bool is_valid_collapse(const Ring& ring, const CollapseCandidate& c);
bool is_valid_collapse_global(const Polygon& poly, int ring_index, const CollapseCandidate& c);

#endif