#include "simplifier.h"
#include "geometry.h"
#include "polygon.h"

#include <cmath>
#include <queue>
#include <vector>
#include <algorithm>
#include <set>

namespace {

/**
 * Axis-Aligned Bounding Box (AABB) for quick intersection pruning.
 */
struct Rect {
    double min_x, min_y, max_x, max_y;

    /**
     * Checks if this rectangle intersects another, including an EPS buffer.
     */
    bool intersects(const Rect& other) const {
        return !(max_x < other.min_x - EPS || min_x > other.max_x + EPS ||
                 max_y < other.min_y - EPS || min_y > other.max_y + EPS);
    }

    /**
     * Creates an AABB from a line segment with an EPS buffer.
     */
    static Rect from_segment(const Point& a, const Point& b) {
        return {std::min(a.x, b.x) - EPS, std::min(a.y, b.y) - EPS,
                std::max(a.x, b.x) + EPS, std::max(a.y, b.y) + EPS};
    }
};

/**
 * A uniform grid spatial index to accelerate segment intersection queries.
 * Divides the polygon's bounding box into a grid and maps segments to cells.
 */
class SpatialIndex {
public:
    /**
     * Builds the grid for all segments in the polygon.
     */
    SpatialIndex(const Polygon& poly, int grid_size = 100) : grid_size_(grid_size) {
        if (poly.rings.empty()) return;

        // Calculate global bounding box
        bool first = true;
        for (const auto& ring : poly.rings) {
            for (const auto& p : ring.vertices) {
                if (first) {
                    bbox_ = {p.x, p.y, p.x, p.y};
                    first = false;
                } else {
                    bbox_.min_x = std::min(bbox_.min_x, p.x);
                    bbox_.min_y = std::min(bbox_.min_y, p.y);
                    bbox_.max_x = std::max(bbox_.max_x, p.x);
                    bbox_.max_y = std::max(bbox_.max_y, p.y);
                }
            }
        }

        // Determine cell dimensions
        cell_w_ = (bbox_.max_x - bbox_.min_x) / grid_size_;
        cell_h_ = (bbox_.max_y - bbox_.min_y) / grid_size_;

        // Avoid division by zero for point or line-like polygons
        if (cell_w_ < EPS) cell_w_ = 1.0;
        if (cell_h_ < EPS) cell_h_ = 1.0;

        grid_.resize(grid_size_ * grid_size_);

        // Populate grid with segments
        for (int r = 0; r < static_cast<int>(poly.rings.size()); ++r) {
            const auto& ring = poly.rings[r];
            for (int i = 0; i < static_cast<int>(ring.vertices.size()); ++i) {
                add_segment(r, i, ring.vertices[i], ring.vertices[(i + 1) % ring.vertices.size()]);
            }
        }
    }

    /**
     * Maps a single segment to all cells it touches.
     */
    void add_segment(int ring_idx, int seg_idx, const Point& a, const Point& b) {
        Rect seg_bbox = Rect::from_segment(a, b);
        int x1 = get_x_cell(seg_bbox.min_x);
        int y1 = get_y_cell(seg_bbox.min_y);
        int x2 = get_x_cell(seg_bbox.max_x);
        int y2 = get_y_cell(seg_bbox.max_y);

        for (int x = x1; x <= x2; ++x) {
            for (int y = y1; y <= y2; ++y) {
                grid_[y * grid_size_ + x].push_back({ring_idx, seg_idx});
            }
        }
    }

    /**
     * Reference to a specific segment within the Polygon.
     */
    struct SegmentRef {
        int ring_idx;
        int seg_idx;
        bool operator<(const SegmentRef& other) const {
            if (ring_idx != other.ring_idx) return ring_idx < other.ring_idx;
            return seg_idx < other.seg_idx;
        }
    };

    /**
     * Returns a set of unique segments that overlap with the given rectangle.
     */
    std::set<SegmentRef> query(const Rect& rect) const {
        std::set<SegmentRef> results;
        int x1 = get_x_cell(rect.min_x);
        int y1 = get_y_cell(rect.min_y);
        int x2 = get_x_cell(rect.max_x);
        int y2 = get_y_cell(rect.max_y);

        for (int x = x1; x <= x2; ++x) {
            for (int y = y1; y <= y2; ++y) {
                const auto& cell = grid_[y * grid_size_ + x];
                for (const auto& ref : cell) {
                    results.insert(ref);
                }
            }
        }
        return results;
    }

private:
    int get_x_cell(double x) const {
        int c = static_cast<int>((x - bbox_.min_x) / cell_w_);
        return std::max(0, std::min(grid_size_ - 1, c));
    }
    int get_y_cell(double y) const {
        int c = static_cast<int>((y - bbox_.min_y) / cell_h_);
        return std::max(0, std::min(grid_size_ - 1, c));
    }

    int grid_size_;
    Rect bbox_{0,0,0,0};
    double cell_w_{1}, cell_h_{1};
    std::vector<std::vector<SegmentRef>> grid_;
};

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
    if (!c.valid || n < 3) return ring;
    Ring out;
    out.ring_id = ring.ring_id;
    const int A = c.start_index % n;
    const int B = (A + 1) % n;
    const int C = (A + 2) % n;
    out.vertices.reserve(n - 1);
    for (int i = 0; i < n; ++i) {
        if (i == B || i == C) continue;
        out.vertices.push_back(ring.vertices[i]);
        if (i == A) out.vertices.push_back(c.e);
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
        if (((a.y > p.y) != (b.y > p.y)) &&
            (p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x)) {
            inside = !inside;
        }
    }
    return inside;
}

/**
 * Checks whether an inner ring lies inside an outer ring using a sample point.
 */
bool ring_inside_ring(const Ring& inner, const Ring& outer) {
    if (inner.vertices.empty()) return false;
    return point_in_ring(inner.vertices[0], outer);
}

} // namespace

/**
 * Calculates the signed area of triangle (a,b,c).
 */
double triangle_signed_area(const Point& a, const Point& b, const Point& c) {
    return 0.5 * ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}

/**
 * Calculates the absolute area of triangle (a,b,c).
 */
double triangle_area_abs(const Point& a, const Point& b, const Point& c) {
    return std::fabs(triangle_signed_area(a, b, c));
}

/**
 * Calculates the local areal displacement for a segment collapse A-B-C-D -> A-E-D.
 * This is the sum of the areas of triangles formed by the removed segments and E.
 */
double exact_local_displacement(const Point& A, const Point& B, const Point& C, const Point& D, const Point& E) {
    return triangle_area_abs(A, B, E) + triangle_area_abs(B, C, E) + triangle_area_abs(C, D, E);
}

/**
 * Computes the optimal candidate point E that preserves the area of segment A-B-C-D.
 * Evaluates candidates on line AB, CD, BC, and the midpoint projection to find the 
 * one with minimum geometric distortion (displacement).
 */
CollapseCandidate compute_collapse_candidate(const Ring& ring, int start_index) {
    CollapseCandidate candidate{};
    const int n = static_cast<int>(ring.vertices.size());
    if (n < 3) return candidate;

    const Point& A = ring.vertices[start_index % n];
    const Point& B = ring.vertices[(start_index + 1) % n];
    const Point& C = ring.vertices[(start_index + 2) % n];
    const Point& D = ring.vertices[(start_index + 3) % n];

    // Total signed area of the polygon region formed by A-B-C-D
    const double area_abcd = triangle_signed_area(A, B, C) + triangle_signed_area(A, C, D);
    const double area_abd = triangle_signed_area(A, B, D);
    const double area_acd = triangle_signed_area(A, C, D);
    const double dx = D.x - A.x;
    const double dy = D.y - A.y;
    const double len_sq = dx * dx + dy * dy;

    if (len_sq < EPS * EPS) return candidate;

    bool found = false;
    Point best_E{};
    double best_disp = 0.0;

    // lambda to update the best candidate based on minimum displacement
    auto update_best = [&](const Point& E) {
        if (same_point(E, A) || same_point(E, D)) return;
        double disp = exact_local_displacement(A, B, C, D, E);
        if (!found || disp < best_disp) {
            best_E = E;
            best_disp = disp;
            found = true;
        }
    };

    // Candidate 1: E on line AB
    if (std::fabs(area_abd) >= EPS) update_best({A.x + (area_abcd / area_abd) * (B.x - A.x), A.y + (area_abcd / area_abd) * (B.y - A.y)});
    
    // Candidate 2: E on line CD
    if (std::fabs(area_acd) >= EPS) update_best({D.x + (area_abcd / area_acd) * (C.x - D.x), D.y + (area_abcd / area_acd) * (C.y - D.y)});
    
    // Candidate 3: E on line BC
    if (std::fabs(area_acd - area_abd) >= EPS) {
        double t = (area_abcd - area_abd) / (area_acd - area_abd);
        update_best({B.x + t * (C.x - B.x), B.y + t * (C.y - B.y)});
    }
    
    // Candidate 4: Perpendicular projection of midpoint BC onto area-preserving line
    Point M{(B.x + C.x) * 0.5, (B.y + C.y) * 0.5};
    double h = 2.0 * (area_abcd - triangle_signed_area(A, M, D)) / std::sqrt(len_sq);
    update_best({M.x - h * (dy / std::sqrt(len_sq)), M.y + h * (dx / std::sqrt(len_sq))});

    if (!found) return candidate;
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
    if (!c.valid) return;
    const int n = static_cast<int>(ring.vertices.size());
    if (n < 3) return;
    const int A = c.start_index % n;
    const int B = (A + 1) % n;
    const int C = (A + 2) % n;
    
    std::vector<Point> new_vertices;
    new_vertices.reserve(n - 1);
    for (int i = 0; i < n; ++i) {
        if (i == B || i == C) continue;
        new_vertices.push_back(ring.vertices[i]);
        if (i == A) new_vertices.push_back(c.e);
    }
    ring.vertices = std::move(new_vertices);
}

/**
 * Global validity check for a collapse candidate.
 * Ensures no intersections with any other parts of the polygon and preserves containment.
 */
bool is_valid_collapse(const Polygon& poly, int ring_idx, const CollapseCandidate& c, const SpatialIndex& index) {
    const Ring& ring = poly.rings[ring_idx];
    const int n = static_cast<int>(ring.vertices.size());
    const int A = c.start_index % n;
    const int B = (A + 1) % n;
    const int C = (A + 2) % n;
    const int D = (A + 3) % n;
    const Point a = ring.vertices[A];
    const Point d = ring.vertices[D];
    const Point e = c.e;

    // Bounding boxes for the two new segments AE and ED
    Rect bb_ae = Rect::from_segment(a, e);
    Rect bb_ed = Rect::from_segment(e, d);
    Rect bb_total = {std::min(bb_ae.min_x, bb_ed.min_x), std::min(bb_ae.min_y, bb_ed.min_y),
                     std::max(bb_ae.max_x, bb_ed.max_x), std::max(bb_ae.max_y, bb_ed.max_y)};

    // Query spatial index for potential intersecting segments
    for (const auto& ref : index.query(bb_total)) {
        const Ring& other_ring = poly.rings[ref.ring_idx];
        const int m = static_cast<int>(other_ring.vertices.size());
        
        // Skip current and adjacent segments in the same ring
        if (ref.ring_idx == ring_idx) {
            int i = ref.seg_idx;
            if (i == (A - 1 + n) % n || i == A || i == B || i == C || i == D) continue;
        }

        const Point& p = other_ring.vertices[ref.seg_idx];
        const Point& q = other_ring.vertices[(ref.seg_idx + 1) % m];
        
        Rect bb_pq = Rect::from_segment(p, q);
        if (bb_ae.intersects(bb_pq)) {
            if (segments_intersect(a, e, p, q)) {
                // Allow touching at endpoints, reject proper internal intersections
                if (!same_point(a, p) && !same_point(a, q) && !same_point(e, p) && !same_point(e, q)) return false;
            }
        }
        if (bb_ed.intersects(bb_pq)) {
            if (segments_intersect(e, d, p, q)) {
                if (!same_point(e, p) && !same_point(e, q) && !same_point(d, p) && !same_point(d, q)) return false;
            }
        }
    }

    // Preserve containment relationships (e.g., holes must stay inside exterior ring)
    Ring changed_ring = make_collapsed_ring(ring, c);
    if (ring_idx == 0) {
        // If exterior ring changed, ensure all holes are still inside
        for (int r = 1; r < static_cast<int>(poly.rings.size()); ++r) {
            if (!ring_inside_ring(poly.rings[r], changed_ring)) return false;
        }
    } else {
        // If a hole changed, ensure it is still inside the exterior ring
        if (!ring_inside_ring(changed_ring, poly.rings[0])) return false;
    }

    return true;
}

/**
 * Main simplification loop using the APSC algorithm.
 * Selects candidates greedily from a priority queue based on minimum displacement.
 */
Polygon simplify_polygon(const Polygon& input, int target_vertices, double& areal_displacement) {
    areal_displacement = 0.0;
    Polygon output = input;
    std::vector<unsigned> ring_versions(output.rings.size(), 0);
    std::priority_queue<PQCandidate, std::vector<PQCandidate>, CandidateCompare> pq;

    // Lambda to queue all possible segment collapses for a given ring
    auto push_ring_candidates = [&](int ring_idx) {
        Ring& ring = output.rings[ring_idx];
        if (ring.vertices.size() < 3) return;
        for (int i = 0; i < static_cast<int>(ring.vertices.size()); ++i) {
            CollapseCandidate c = compute_collapse_candidate(ring, i);
            if (c.valid) pq.push(PQCandidate{c.areal_displacement, ring_idx, i, ring_versions[ring_idx]});
        }
    };

    for (int r = 0; r < static_cast<int>(output.rings.size()); ++r) push_ring_candidates(r);
    
    SpatialIndex index(output);
    int total_v = total_vertex_count(output);

    while (total_v > target_vertices && !pq.empty()) {
        PQCandidate top = pq.top();
        pq.pop();
        
        // Skip stale candidates (from modified rings)
        if (top.ring_idx < 0 || top.ring_idx >= static_cast<int>(output.rings.size()) || top.ring_version != ring_versions[top.ring_idx]) continue;
        
        Ring& ring = output.rings[top.ring_idx];
        if (ring.vertices.size() < 3) continue;
        
        CollapseCandidate c = compute_collapse_candidate(ring, top.start_index);
        
        // Apply collapse if topology is preserved
        if (c.valid && is_valid_collapse(output, top.ring_idx, c, index)) {
            apply_collapse(ring, c);
            areal_displacement += c.areal_displacement;
            ring_versions[top.ring_idx]++;
            total_v--;
            
            // Periodically rebuild spatial index to maintain accuracy
            if (total_v % 50 == 0) index = SpatialIndex(output);
            
            // Re-queue new candidates for the affected area
            const int new_n = static_cast<int>(ring.vertices.size());
            for (int offset = -3; offset <= 3; ++offset) {
                int idx = (top.start_index + offset + new_n) % new_n;
                CollapseCandidate new_c = compute_collapse_candidate(ring, idx);
                if (new_c.valid) pq.push(PQCandidate{new_c.areal_displacement, top.ring_idx, idx, ring_versions[top.ring_idx]});
            }
        }
    }
    return output;
}
