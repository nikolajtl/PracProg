#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <ranges>

// Task 1
int main(int argc, char* argv[]) {
    std::vector<double> indices; // Defines lists and leaves them empty for appending
    for (int i=0; i<argc; i++) {
        std::string arg = argv[i];
        if (arg=="-n" && i+1<argc) {
            indices.push_back(i+1);
        }
    }

    for (int i: indices) {
        std::string dummy = argv[i];
        double number = std::stod(dummy); // "stod" means "string to double"
        std::cout
        << "Argument #" << i << ": " << number
        << ". sin(" << number << ") = " << std::sin(number)
        << ", cos(" << number << ") = " << std::cos(number)
        << std::endl;
    }
}