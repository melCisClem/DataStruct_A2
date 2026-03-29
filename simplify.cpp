#include "geometry.h"
#include "polygon.h"
#include "simplifier.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: ./simplify <input_file.csv> <target_vertices>\n";
        return 1;
    }

    try {
        const std::string input_file = argv[1];
        const int target_vertices = std::stoi(argv[2]);

        Polygon input_polygon = read_polygon_csv(input_file);
        const double input_area = total_signed_area(input_polygon);

        double areal_displacement = 0.0;
        Polygon output_polygon = simplify_polygon(input_polygon, target_vertices, areal_displacement);
        const double output_area = total_signed_area(output_polygon);

        print_polygon_csv(output_polygon);

        std::cout << std::scientific << std::setprecision(6);
        std::cout << "Total signed area in input: " << input_area << "\n";
        std::cout << "Total signed area in output: " << output_area << "\n";
        std::cout << "Total areal displacement: " << areal_displacement << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}