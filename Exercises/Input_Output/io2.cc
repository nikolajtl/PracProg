#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <ranges>

// Task 2
int main() {
    double input;
    while (std::cin >> input) { // Doesn't terminate automatically
        std::cout
        << "Input: " << input
        << ". sin(" << input << ") = " << std::sin(input)
        << ", cos(" << input << ") = " << std::cos(input)
        << std::endl;
    }
    
    return 0;
}