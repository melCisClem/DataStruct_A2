#include "polygon.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

/**
 * Reads a polygon from a CSV file.
 * The CSV is expected to have a header: ring_id,vertex_id,x,y.
 * Uses a map to group vertices by ring_id and ensures they are stored in the original order.
 */
Polygon read_polygon_csv(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin) {
        throw std::runtime_error("Failed to open input file: " + filename);
    }

    std::string line;
    // Skip header line
    if (!std::getline(fin, line)) {
        throw std::runtime_error("Input file is empty.");
    }

    // Temporary map to collect rings by their ID
    std::map<int, Ring> ring_map;

    while (std::getline(fin, line)) {
        if (line.empty()) {
            continue;
        }

        std::stringstream ss(line);
        std::string token;

        int ring_id = 0;
        double x = 0.0;
        double y = 0.0;

        // Parse: ring_id, vertex_id, x, y
        std::getline(ss, token, ',');
        ring_id = std::stoi(token);

        std::getline(ss, token, ','); // vertex_id (unused, assumed sequential)
        
        std::getline(ss, token, ',');
        x = std::stod(token);

        std::getline(ss, token, ',');
        y = std::stod(token);

        // Initialize new ring if not encountered before
        if (ring_map.find(ring_id) == ring_map.end()) {
            Ring ring;
            ring.ring_id = ring_id;
            ring_map[ring_id] = ring;
        }

        ring_map[ring_id].vertices.push_back({x, y});
    }

    // Convert map to Polygon structure (vector of rings)
    Polygon polygon;
    for (const auto& [id, ring] : ring_map) {
        polygon.rings.push_back(ring);
    }

    return polygon;
}

/**
 * Prints the polygon vertices in CSV format to standard output.
 * Vertex IDs are regenerated sequentially starting from 0 for each ring.
 */
void print_polygon_csv(const Polygon& polygon) {
    std::cout << "ring_id,vertex_id,x,y\n";

    for (const Ring& ring : polygon.rings) {
        for (std::size_t i = 0; i < ring.vertices.size(); ++i) {
            std::cout << ring.ring_id << ","
                      << i << ","
                      << ring.vertices[i].x << ","
                      << ring.vertices[i].y << "\n";
        }
    }
}

/**
 * Returns the total count of vertices across all rings in the polygon.
 */
int total_vertex_count(const Polygon& polygon) {
    int total = 0;
    for (const Ring& ring : polygon.rings) {
        total += static_cast<int>(ring.vertices.size());
    }
    return total;
}
