#pragma once

#include <vector>
#include <functional>
#include <utility> // for std::pair
#include "matrix.h" // assumes pp::vector is defined here

namespace ode {

std::pair<pp::vector, pp::vector> rkstepXY(
    std::function<pp::vector(double, pp::vector)> f,
    double x,
    pp::vector y,
    double h
);

std::pair<std::vector<double>, std::vector<pp::vector>> driver(
    std::function<pp::vector(double, pp::vector)> f,
    std::pair<double, double> interval,
    pp::vector y_a,
    double h = 0.1,
    double acc = 0.01,
    double eps = 0.01
);

} // namespace ode
