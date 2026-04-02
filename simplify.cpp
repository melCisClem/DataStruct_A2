#include "geometry.h"
#include "polygon.h"
#include "simplifier.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

/**
 * Entry point for the polygon simplification tool.
 * Command-line arguments:
 *   1. Path to the input CSV file.
 *   2. Target maximum number of vertices for the simplified polygon.
 */
int main(int argc, char* argv[]) {
    // Validate command-line arguments
    if (argc != 3) {
        std::cerr << "Usage: ./simplify <input_file.csv> <target_vertices>\n";
        return 1;
    }

    try {
        const std::string input_file = argv[1];
        const int target_vertices = std::stoi(argv[2]);

        // Load the polygon from the CSV file
        Polygon input_polygon = read_polygon_csv(input_file);
        
        // Calculate the initial signed area for verification
        const double input_area = total_signed_area(input_polygon);

        // Run the Area-Preserving Segment Collapse (APSC) simplification algorithm
        double areal_displacement = 0.0;
        Polygon output_polygon = simplify_polygon(input_polygon, target_vertices, areal_displacement);
        
        // Calculate final area to ensure area preservation was successful
        const double output_area = total_signed_area(output_polygon);

        // Output the simplified polygon in the specified CSV format
        print_polygon_csv(output_polygon);

        // Output summary information in scientific notation with 6 decimal places
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "Total signed area in input: " << input_area << "\n";
        std::cout << "Total signed area in output: " << output_area << "\n";
        std::cout << "Total areal displacement: " << areal_displacement << "\n";
    }
    catch (const std::exception& e) {
        // Handle file I/O errors or invalid input formats
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
