#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional> // For std::function
#include <algorithm> // For std::min
#include <utility> // For std::pair
#include "matrix.h"

double pi = 3.14159265358979323846;


std::pair< pp::vector, pp::vector > rkstepXY (
    std::function<pp::vector(double, pp::vector)> f, // Driving function f(x, y)
    double x, pp::vector y, double h
) {
    double alpha = 1.0/3.0, beta = 2.0/3.0; // Heun's third order
    pp::vector k0 = f(x, y);
    pp::vector k1 = f(x+alpha*h, y + alpha*k0*h);
    pp::vector k2 = f(x+beta*h, y + (beta-beta*beta/(2*alpha))*k0*h + beta*beta/(2*alpha)*k1*h);
    pp::vector y_next = y + (1 - 1/(2*alpha) - 1/(3*beta*beta) + 1/(3*alpha*beta))*k0*h
                        + (1/(2*alpha) - 1/(3*alpha*beta) + 1/(3*beta*beta))*k1*h;
    pp::vector y_lower = y + (1 - 1/(2*alpha))*k0*h + 1/(2*alpha)*k1*h; // Lower order for error estime 
    pp::vector δy = ( y_next - y_lower )*h; // SIGNED vector of errors
    return {y_next, δy};
};



std::pair< std::vector<double>, std::vector<pp::vector> > driver ( // Want to return xs, ys
    std::function<pp::vector(double, pp::vector)> f, // Driving function f(x, y)
    std::pair<double, double> interval, // (a, b)
    pp::vector y_a, // y(a)=y_a
    double h = 0.1, // Initial step size
    double acc = 0.01, // Absolute accuracy
    double eps = 0.01 // Relative accuracy
) {
    std::vector<double> xs;
    std::vector<pp::vector> ys;
    double a = interval.first;
    double b = interval.second;
    xs.push_back(a);
    ys.push_back(y_a);

    double x = a;
    pp::vector y=y_a;
    do {
        if (x >= b) return {xs, ys}; // Compiler will convert to std::pair from function return type
        if (x+h >= b) x=b; // Should fit with the interval endpoint
        auto [yh, δy] = rkstepXY(f, x, y, h);
        double tol = (acc+eps*yh.norm() * std::sqrt(h/(b-a)));
        double err = δy.norm(); // Norm ensures positive value
        if (err<=tol) {
            x += h; y=yh;
            xs.push_back(x); ys.push_back(y);
        } else { h *= std::min(std::pow(tol/err, 0.25)*0.95, 2.0); }; // Max step size 2
    } while (true); // Careful
    // return {xs, ys};
};


// Exponential decay
pp::vector f1 (double x, pp::vector y) {
    return -1*y;
};

// Harmonic oscillator
pp::vector f2(double x, pp::vector y) {
    // y[0] = y, y[1] = v = y'
    pp::vector result(2);
    result[0] = y[1];      // y' = v
    result[1] = -y[0];     // v' = -y
    return result;
};

// Harmonic oscillator with friction
pp::vector f3(double x, pp::vector y) {
    double eta = 0.1; // Friction coefficient, probably want a small one
    pp::vector result(2);
    result[0] = y[1]; // y_1' = y_2
    result[1] = -y[0] -eta*y[1]; // y_2' = - y_1 - y_2
    return result;
};

// Part B: Relativistic precession
pp::vector f(double x, pp::vector y) {
    double eps = 0.01; // Relativistic correction
    pp::vector result(2);
    result[0] = y[1];
    result[1] = 1 - y[0] + eps*y[0]*y[0];
    return result;
};



int main() {
    std::pair<double, double> interval = {0.0, 8*2*pi}; // Interval
    pp::vector y_a({1, -0.9}); // Initial conditions
    auto [xs, ys] = driver(f, interval, y_a, 0.01, 1e-9, 1e-9);

    std::ofstream funcfile("funcvals.txt");
    if (!funcfile) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    };

    // for (int i=0; i < xs.size(); i++) {
    //     funcfile << xs[i] << " " << ys[i][0] << std::endl;
    // };

    // Transform the coordinates to a polar plot
    for (int i=0; i < xs.size(); i++) {
        funcfile << (1/ys[i][0])*std::cos(xs[i]) << " " << (1/ys[i][0])*std::sin(xs[i]) << std::endl;
    };

    funcfile.close();
    std::cout << "Diff. eqn. solution written to funcvals.txt" << std::endl;

    return 0;
};