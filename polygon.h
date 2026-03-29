#ifndef POLYGON_H
#define POLYGON_H

#include <string>
#include <vector>

struct Point {
    double x{};
    double y{};
};

struct Ring {
    int ring_id{};
    std::vector<Point> vertices;
};

struct Polygon {
    std::vector<Ring> rings;
};

Polygon read_polygon_csv(const std::string& filename);
void print_polygon_csv(const Polygon& polygon);
int total_vertex_count(const Polygon& polygon);

#endif