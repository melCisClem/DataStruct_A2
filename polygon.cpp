#include <iomanip>
#include "polygon.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

/**
 * Reads a polygon from a CSV file.
 * The CSV is expected to have a header: ring_id,vertex_id,x,y.
 * 
 * Vertices are grouped into rings based on their ring_id. The ring_id is assumed
 * to be an integer, and vertices within each ring are assumed to be listed in 
 * sequential order in the input file.
 */
Polygon read_polygon_csv(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin) {
        throw std::runtime_error("Failed to open input file: " + filename);
    }

    std::string line;
    // Skip the header line (ring_id,vertex_id,x,y)
    if (!std::getline(fin, line)) {
        throw std::runtime_error("Input file is empty.");
    }

    // Temporary map to collect rings by their ID. std::map keeps IDs sorted.
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

        // Extract ring_id
        if (!std::getline(ss, token, ',')) continue;
        ring_id = std::stoi(token);

        // vertex_id is skipped as we regenerate it sequentially on output
        if (!std::getline(ss, token, ',')) continue;
        
        // Extract x coordinate
        if (!std::getline(ss, token, ',')) continue;
        x = std::stod(token);

        // Extract y coordinate
        if (!std::getline(ss, token, ',')) continue;
        y = std::stod(token);

        // Initialize new ring if this ID hasn't been seen yet
        if (ring_map.find(ring_id) == ring_map.end()) {
            Ring ring;
            ring.ring_id = ring_id;
            ring_map[ring_id] = ring;
        }

        // Append vertex to the appropriate ring
        ring_map[ring_id].vertices.push_back({x, y});
    }

    // Transfer rings from the sorted map to the Polygon's vector
    Polygon polygon;
    for (const auto& [id, ring] : ring_map) {
        polygon.rings.push_back(ring);
    }

    return polygon;
}

/**
 * Prints the polygon vertices in CSV format to standard output.
 * 
 * Vertex coordinates are printed with high precision (10 decimal places) 
 * to ensure consistency with calculated area and preservation metrics.
 * Vertex IDs are regenerated sequentially (0, 1, 2, ...) within each ring.
 */
void print_polygon_csv(const Polygon& polygon) {
    std::cout << "ring_id,vertex_id,x,y\n";
    // Set higher precision for coordinate output to ensure consistency with summary area.
    // 10 decimal places should be sufficient for most cases.
    std::cout << std::setprecision(10);

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
 * Used to check against the target vertex count during simplification.
 */
int total_vertex_count(const Polygon& polygon) {
    int total = 0;
    for (const Ring& ring : polygon.rings) {
        total += static_cast<int>(ring.vertices.size());
    }
    return total;
}
