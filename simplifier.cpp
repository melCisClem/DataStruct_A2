#include "simplifier.h"
#include "geometry.h"
#include "polygon.h"

#include <cmath>
#include <iostream>
#include <queue>
#include <vector>

namespace {

constexpr bool DEBUG_LOG = true;

struct PQCandidate {
    double displacement{0.0};
    int ring_idx{-1};
    int start_index{-1};
    unsigned ring_version{0};
};

struct CandidateCompare {
    bool operator()(const PQCandidate& a, const PQCandidate& b) const {
        return a.displacement > b.displacement; // min-heap behavior
    }
};

Ring make_collapsed_ring(const Ring& ring, const CollapseCandidate& c) {
    Ring out;
    out.ring_id = ring.ring_id;

    const int n = static_cast<int>(ring.vertices.size());
    if (!c.valid || n < 4) {
        return ring;
    }

    const int A = c.start_index % n;
    const int B = (A + 1) % n;
    const int C = (A + 2) % n;

    out.vertices.reserve(n - 1);

    for (int i = 0; i < n; ++i) {
        if (i == B || i == C) {
            continue;
        }

        if (i == A) {
            out.vertices.push_back(ring.vertices[A]);
            out.vertices.push_back(c.e);
        } else {
            out.vertices.push_back(ring.vertices[i]);
        }
    }

    return out;
}

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

bool ring_inside_ring(const Ring& inner, const Ring& outer) {
    for (const Point& p : inner.vertices) {
        if (!point_in_ring(p, outer)) {
            return false;
        }
    }
    return true;
}

double exact_local_displacement(const Point& A,
                                const Point& B,
                                const Point& C,
                                const Point& D,
                                const Point& E) {
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

    const double target_signed_area =
        triangle_signed_area(A, B, C) + triangle_signed_area(A, C, D);

    const double dx = D.x - A.x;
    const double dy = D.y - A.y;
    const double len = std::sqrt(dx * dx + dy * dy);

    if (len < EPS) {
        return candidate;
    }

    const double nx = -dy / len;
    const double ny = dx / len;

    // triangle_signed_area(A, E, D) = -0.5 * signed_height * len
    const double signed_height = (-2.0 * target_signed_area) / len;

    Point E{
        A.x + signed_height * nx,
        A.y + signed_height * ny
    };

    candidate.valid = true;
    candidate.e = E;
    candidate.start_index = start_index;
    candidate.areal_displacement = exact_local_displacement(A, B, C, D, E);

    if (DEBUG_LOG) {
        std::cerr << "[DEBUG] ring " << ring.ring_id
                  << " target_signed_area=" << target_signed_area
                  << " |AD|=" << len << "\n";
        std::cerr << "[DEBUG] ring " << ring.ring_id
                  << " E=(" << E.x << ", " << E.y << ")"
                  << " displacement=" << candidate.areal_displacement << "\n";
    }

    return candidate;
}

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

        if (i == A) {
            new_vertices.push_back(ring.vertices[A]);
            new_vertices.push_back(c.e);
        } else {
            new_vertices.push_back(ring.vertices[i]);
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
        const int j = (i + 1) % n;

        // Skip local edges being replaced
        if (i == A || i == B || i == C) continue;
        if (j == A || j == B || j == C) continue;

        const Point p = ring.vertices[i];
        const Point q = ring.vertices[j];

        // Skip edges sharing endpoints with new edges
        if (same_point(p, a) || same_point(p, d) ||
            same_point(q, a) || same_point(q, d)) {
            continue;
        }

        if (segments_intersect(a, e, p, q)) return false;
        if (segments_intersect(e, d, p, q)) return false;
    }

    return true;
}

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

    // 1. same-ring checks
    if (!is_valid_collapse(ring, c)) {
        return false;
    }

    // 2. cross-ring intersection checks
    for (int r = 0; r < static_cast<int>(poly.rings.size()); ++r) {
        if (r == ring_idx) continue;

        const Ring& other = poly.rings[r];
        const int m = static_cast<int>(other.vertices.size());

        for (int i = 0; i < m; ++i) {
            const int j = (i + 1) % m;
            const Point& p = other.vertices[i];
            const Point& q = other.vertices[j];

            if (segments_intersect(a, e, p, q)) return false;
            if (segments_intersect(e, d, p, q)) return false;
        }
    }

    // 3. containment checks using the UPDATED ring
    Ring changed_ring = make_collapsed_ring(ring, c);

    if (ring_idx == 0) {
        // Outer ring changed: every hole must still lie inside it
        for (int r = 1; r < static_cast<int>(poly.rings.size()); ++r) {
            if (!ring_inside_ring(poly.rings[r], changed_ring)) {
                return false;
            }
        }
    } else {
        // Hole changed: changed hole must still lie inside outer ring
        if (!ring_inside_ring(changed_ring, poly.rings[0])) {
            return false;
        }
    }

    return true;
}

Polygon simplify_polygon(const Polygon& input, int target_vertices, double& areal_displacement) {
    areal_displacement = 0.0;

    Polygon output = input;
    std::vector<unsigned> ring_versions(output.rings.size(), 0);

    if (DEBUG_LOG) {
        std::cerr << "[DEBUG] Target vertices: " << target_vertices << "\n";
        std::cerr << "[DEBUG] Total input vertices: " << total_vertex_count(input) << "\n";
    }

    std::priority_queue<PQCandidate, std::vector<PQCandidate>, CandidateCompare> pq;

    auto push_ring_candidates = [&](int ring_idx) {
        Ring& ring = output.rings[ring_idx];
        if (ring.vertices.size() < 4) {
            return;
        }

        for (int i = 0; i < static_cast<int>(ring.vertices.size()); ++i) {
            CollapseCandidate c = compute_collapse_candidate(ring, i);
            if (!c.valid) continue;

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
            continue; // stale candidate
        }

        Ring& ring = output.rings[top.ring_idx];
        if (ring.vertices.size() < 4) {
            continue;
        }

        CollapseCandidate c = compute_collapse_candidate(ring, top.start_index);
        if (!c.valid) {
            continue;
        }

        if (!is_valid_collapse_global(output, top.ring_idx, c)) {
            continue;
        }

        if (DEBUG_LOG) {
            std::cerr << "[DEBUG] Applying best collapse on ring "
                      << ring.ring_id
                      << " with displacement = " << c.areal_displacement
                      << "\n";
        }

        apply_collapse(ring, c);
        areal_displacement += c.areal_displacement;
        ring_versions[top.ring_idx]++;

        if (DEBUG_LOG) {
            std::cerr << "[DEBUG] Total vertices after collapse: "
                      << total_vertex_count(output) << "\n";
        }

        // local-ish update: rebuild candidates for affected ring only
        push_ring_candidates(top.ring_idx);
    }

    if (DEBUG_LOG && total_vertex_count(output) > target_vertices) {
        std::cerr << "[DEBUG] No valid collapse found. Stopping.\n";
    }

    return output;
}