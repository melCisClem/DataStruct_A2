#include "simplifier.h"
#include "geometry.h"
#include "polygon.h"

#include <cmath>
#include <queue>
#include <vector>

namespace {

/**
 * Representation of a candidate collapse in the Priority Queue.
 * Includes versioning to handle stale candidates after a ring is modified.
 */
struct PQCandidate {
    double displacement{0.0};
    int ring_idx{-1};
    int start_index{-1};
    unsigned ring_version{0};
};

/**
 * Comparator for the Min-Heap priority queue to always pick the collapse
 * with the smallest areal displacement.
 */
struct CandidateCompare {
    bool operator()(const PQCandidate& a, const PQCandidate& b) const {
        if (std::fabs(a.displacement - b.displacement) > EPS) {
            return a.displacement > b.displacement;
        }

        if (a.ring_idx != b.ring_idx) {
            return a.ring_idx > b.ring_idx;
        }

        return a.start_index > b.start_index;
    }
};

/**
 * Helper to create a new ring object with a collapse applied.
 * Replaces vertices at (A+1) and (A+2) with point E.
 */
Ring make_collapsed_ring(const Ring& ring, const CollapseCandidate& c) {
    const int n = static_cast<int>(ring.vertices.size());
    if (!c.valid || n < 4) {
        return ring;
    }

    Ring out;
    out.ring_id = ring.ring_id;

    const int A = c.start_index % n;
    const int B = (A + 1) % n;
    const int C = (A + 2) % n;

    out.vertices.reserve(n - 1);

    for (int i = 0; i < n; ++i) {
        if (i == B || i == C) {
            continue;
        }

        out.vertices.push_back(ring.vertices[i]);
        if (i == A) {
            out.vertices.push_back(c.e);
        }
    }

    return out;
}

/**
 * Standard ray-casting algorithm to check if a point lies inside a ring.
 */
bool point_in_ring(const Point& p, const Ring& ring) {
    bool inside = false;
    const int n = static_cast<int>(ring.vertices.size());

    for (int i = 0, j = n - 1; i < n; j = i++) {
        const Point& a = ring.vertices[i];
        const Point& b = ring.vertices[j];

        const bool crosses =
            ((a.y > p.y) != (b.y > p.y)) &&
            (p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x);

        if (crosses) {
            inside = !inside;
        }
    }

    return inside;
}

/**
 * Checks whether an inner ring lies inside an outer ring.
 * A single point is enough here because rings are assumed valid/simple and
 * intersection checks are handled elsewhere.
 */
bool ring_inside_ring(const Ring& inner, const Ring& outer) {
    if (inner.vertices.empty()) {
        return false;
    }
    return point_in_ring(inner.vertices[0], outer);
}

/**
 * Calculates the local areal displacement for a segment collapse A-B-C-D -> A-E-D.
 */
double exact_local_displacement(const Point& A,
                                const Point& B,
                                const Point& C,
                                const Point& D,
                                const Point& E) {
    if (std::fabs(cross(A, B, E)) < EPS) {
        return triangle_area_abs(E, B, C) +
               triangle_area_abs(E, C, D);
    }

    if (std::fabs(cross(C, D, E)) < EPS) {
        return triangle_area_abs(A, B, E) +
               triangle_area_abs(B, C, E);
    }

    return triangle_area_abs(A, B, E) +
           triangle_area_abs(B, C, E) +
           triangle_area_abs(C, D, E);
}

} // namespace

double triangle_signed_area(const Point& a, const Point& b, const Point& c) {
    return 0.5 * ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}

double triangle_area_abs(const Point& a, const Point& b, const Point& c) {
    return std::fabs(triangle_signed_area(a, b, c));
}

/**
 * Computes a candidate point E that preserves the area of the segment A-B-C-D.
 */
CollapseCandidate compute_collapse_candidate(const Ring& ring, int start_index) {
    CollapseCandidate candidate{};
    const int n = static_cast<int>(ring.vertices.size());

    if (n < 4) {
        return candidate;
    }

    const Point& A = ring.vertices[start_index % n];
    const Point& B = ring.vertices[(start_index + 1) % n];
    const Point& C = ring.vertices[(start_index + 2) % n];
    const Point& D = ring.vertices[(start_index + 3) % n];

    const double area_abcd =
        triangle_signed_area(A, B, C) +
        triangle_signed_area(A, C, D);

    const double area_abd = triangle_signed_area(A, B, D);
    const double area_acd = triangle_signed_area(A, C, D);

    const double dx = D.x - A.x;
    const double dy = D.y - A.y;
    const double len_sq = dx * dx + dy * dy;

    if (len_sq < EPS * EPS) {
        return candidate;
    }

    bool found = false;
    Point best_E{};
    double best_disp = 0.0;

    // Candidate 1: E lies on the LINE through AB
    if (std::fabs(area_abd) >= EPS) {
        const double t_ab = area_abcd / area_abd;

        Point E_ab{
            A.x + t_ab * (B.x - A.x),
            A.y + t_ab * (B.y - A.y)
        };

        const double disp_ab = exact_local_displacement(A, B, C, D, E_ab);
        best_E = E_ab;
        best_disp = disp_ab;
        found = true;
    }

    // Candidate 2: E lies on the LINE through CD
    if (std::fabs(area_acd) >= EPS) {
        const double t_cd = area_abcd / area_acd;

        Point E_cd{
            D.x + t_cd * (C.x - D.x),
            D.y + t_cd * (C.y - D.y)
        };

        const double disp_cd = exact_local_displacement(A, B, C, D, E_cd);

        if (!found || disp_cd < best_disp) {
            best_E = E_cd;
            best_disp = disp_cd;
            found = true;
        }
    }

    if (!found) {
        return candidate;
    }

    candidate.valid = true;
    candidate.e = best_E;
    candidate.start_index = start_index;
    candidate.areal_displacement = best_disp;

    return candidate;
}

/**
 * Modifies the ring in-place to apply the chosen collapse.
 */
void apply_collapse(Ring& ring, const CollapseCandidate& c) {
    if (!c.valid) {
        return;
    }

    const int n = static_cast<int>(ring.vertices.size());
    if (n < 4) {
        return;
    }

    const int A = c.start_index % n;
    const int B = (A + 1) % n;
    const int C = (A + 2) % n;

    std::vector<Point> new_vertices;
    new_vertices.reserve(n - 1);

    for (int i = 0; i < n; ++i) {
        if (i == B || i == C) {
            continue;
        }

        new_vertices.push_back(ring.vertices[i]);
        if (i == A) {
            new_vertices.push_back(c.e);
        }
    }

    ring.vertices = std::move(new_vertices);
}

bool is_valid_collapse(const Ring& ring, const CollapseCandidate& c) {
    const int n = static_cast<int>(ring.vertices.size());
    if (n < 4 || !c.valid) {
        return false;
    }

    const int A = c.start_index % n;
    const int B = (A + 1) % n;
    const int C = (A + 2) % n;
    const int D = (A + 3) % n;

    const Point a = ring.vertices[A];
    const Point d = ring.vertices[D];
    const Point e = c.e;

    for (int i = 0; i < n; ++i) {
        // Skip edges adjacent to the collapsed region:
        // (A-1,A), (A,B), (B,C), (C,D), (D,D+1)
        if (i == (A - 1 + n) % n || i == A || i == B || i == C || i == D) {
            continue;
        }

        const int j = (i + 1) % n;
        const Point& p = ring.vertices[i];
        const Point& q = ring.vertices[j];

        if (segments_intersect(a, e, p, q)) {
            return false;
        }
        if (segments_intersect(e, d, p, q)) {
            return false;
        }
    }

    return true;
}

/**
 * Global validity check for a collapse candidate.
 * Ensures no self-intersection, no intersection with other rings,
 * and preserves containment.
 */
bool is_valid_collapse_global(const Polygon& poly, int ring_idx, const CollapseCandidate& c) {
    const Ring& ring = poly.rings[ring_idx];
    const int n = static_cast<int>(ring.vertices.size());

    if (n < 4 || !c.valid) {
        return false;
    }

    const int A = c.start_index % n;
    const int D = (A + 3) % n;

    const Point a = ring.vertices[A];
    const Point d = ring.vertices[D];
    const Point e = c.e;

    if (!is_valid_collapse(ring, c)) {
        return false;
    }

    for (int r = 0; r < static_cast<int>(poly.rings.size()); ++r) {
        if (r == ring_idx) {
            continue;
        }

        const Ring& other = poly.rings[r];
        const int m = static_cast<int>(other.vertices.size());

        for (int i = 0; i < m; ++i) {
            const int j = (i + 1) % m;
            const Point& p = other.vertices[i];
            const Point& q = other.vertices[j];

            if (segments_intersect(a, e, p, q)) {
                return false;
            }
            if (segments_intersect(e, d, p, q)) {
                return false;
            }
        }
    }

    Ring changed_ring = make_collapsed_ring(ring, c);

    if (ring_idx == 0) {
        for (int r = 1; r < static_cast<int>(poly.rings.size()); ++r) {
            if (!ring_inside_ring(poly.rings[r], changed_ring)) {
                return false;
            }
        }
    } else {
        if (!ring_inside_ring(changed_ring, poly.rings[0])) {
            return false;
        }
    }

    return true;
}

/**
 * Main simplification loop.
 */
Polygon simplify_polygon(const Polygon& input, int target_vertices, double& areal_displacement) {
    areal_displacement = 0.0;

    Polygon output = input;
    std::vector<unsigned> ring_versions(output.rings.size(), 0);

    std::priority_queue<PQCandidate, std::vector<PQCandidate>, CandidateCompare> pq;

    auto push_ring_candidates = [&](int ring_idx) {
        Ring& ring = output.rings[ring_idx];
        if (ring.vertices.size() <= 3) {
            return;
        }

        for (int i = 0; i < static_cast<int>(ring.vertices.size()); ++i) {
            CollapseCandidate c = compute_collapse_candidate(ring, i);
            if (!c.valid) {
                continue;
            }

            pq.push(PQCandidate{
                c.areal_displacement,
                ring_idx,
                i,
                ring_versions[ring_idx]
            });
        }
    };

    for (int r = 0; r < static_cast<int>(output.rings.size()); ++r) {
        push_ring_candidates(r);
    }

    while (total_vertex_count(output) > target_vertices && !pq.empty()) {
        PQCandidate top = pq.top();
        pq.pop();

        if (top.ring_idx < 0 || top.ring_idx >= static_cast<int>(output.rings.size())) {
            continue;
        }

        if (top.ring_version != ring_versions[top.ring_idx]) {
            continue;
        }

        Ring& ring = output.rings[top.ring_idx];
        if (ring.vertices.size() <= 3) {
            continue;
        }

        CollapseCandidate c = compute_collapse_candidate(ring, top.start_index);
        if (!c.valid) {
            continue;
        }

        if (!is_valid_collapse(ring, c)) {
            continue;
        }

        if (!is_valid_collapse_global(output, top.ring_idx, c)) {
            continue;
        }

        apply_collapse(ring, c);
        areal_displacement += c.areal_displacement;
        ring_versions[top.ring_idx]++;

        const int new_n = static_cast<int>(ring.vertices.size());
        for (int offset = -3; offset <= 3; ++offset) {
            int idx = (top.start_index + offset + new_n) % new_n;
            CollapseCandidate new_c = compute_collapse_candidate(ring, idx);
            if (!new_c.valid) {
                continue;
            }

            pq.push(PQCandidate{
                new_c.areal_displacement,
                top.ring_idx,
                idx,
                ring_versions[top.ring_idx]
            });
        }
    }

    return output;
}