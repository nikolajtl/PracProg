#include <iostream>
#include <fstream> // For writing out the data
#include <vector>
#include <cmath>
#include <random>
#include <functional> // For std::function
#include <algorithm> // For std::transform
#include "matrix.h"


std::vector<double> linspace(double start, double end, int num_points) {
    std::vector<double> result;
    double step = (end - start) / (num_points - 1);  // Calculate step size

    for (int i = 0; i < num_points; ++i) {
        result.push_back(start + i * step);
    }
    return result;
}


int binsearch(double x, std::vector<double> xs) {
    if (x < xs[0] || x>xs[xs.size()-1]) throw std::runtime_error("Binary search target out list range.");
    int i=0, j=xs.size()-1;
    while (j-i>1) {
        int av = (i+j) / 2;
        if (x>xs[av]) i = av; else j=av; // Takes the lower bound
    };
    return i;
};


double linterp(double z, std::vector<double> xs, std::vector<double> ys) { // Produces the value val=f_linspline(z) given the lists of (x_i, y_i) points
    int i=binsearch(z, xs);
    double dx = xs[i+1] - xs[i];
    double dy = ys[i+1] - ys[i];
    return ys[i]+dy/dx*(z-xs[i]);
};


double lintegral(double z, std::vector<double> xs, std::vector<double> ys) {
    int i=binsearch(z, xs);
    double integral = 0;
    for (int j=0; j<i; j++) {
        double dx = xs[j+1] - xs[j];
        double dy = ys[j+1] - ys[j];
        double p = dy/dx;
        integral += ys[j]*dx + p*dx*dx / 2;
    };
    double dx = z - xs[i];
    double dy = linterp(z, xs, ys) - ys[i];
    double p = dy/dx;
    integral += ys[i]*dx + p*dx*dx / 2;
    return integral;
};


double spline(double x, std::vector<double> xs, std::vector<double> ys) {
    double y = linterp(x, xs, ys);
    return y;
};


int main() {
    // Dummy x-vals
    std::vector<double> xs = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    // Make cos(x) points
    std::vector<double> ys(xs.size());
    std::transform(xs.begin(), xs.end(), ys.begin(), [](double n) { return std::cos(n); });


    // Write xs and ys to data.txt
    std::ofstream datfile("data.txt");
    if (!datfile) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    for (size_t i = 0; i < xs.size(); ++i) {
        datfile << xs[i] << " " << ys[i] << std::endl;
    }

    datfile.close();
    std::cout << "Example points written to data.txt" << std::endl;


    // Write a spline to splinedat.txt
    std::ofstream splinefile("splinedata.txt");
    if (!splinefile) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // Make the spline x-space, y-space and y-space of the integral
    std::vector<double> xs_spline = linspace(0, 10, 100);

    std::vector<double> ys_spline(xs_spline.size());
    std::vector<double> ys_int(xs_spline.size());
    std::transform(xs_spline.begin(), xs_spline.end(), ys_spline.begin(),
                   [&](double x) { return spline(x, xs, ys); });
    std::transform(xs_spline.begin(), xs_spline.end(), ys_int.begin(),
                   [&](double x) { return lintegral(x, xs, ys); });
    
    for (size_t i = 0; i < xs_spline.size(); ++i) {
        splinefile << xs_spline[i] << " " << ys_spline[i] << std::endl;
    }

    splinefile.close();
    std::cout << "Spline data written to splinedata.txt" << std::endl;

    // Write the integral into a file
    std::ofstream intfile("integraldata.txt");
    if (!intfile) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    for (int i=0; i < xs_spline.size(); i++) {
        intfile << xs_spline[i] << " " << ys_int[i] << std::endl;
    }

    intfile.close();
    std::cout << "Spline integral data written to integraldata.txt" << std::endl;

    return 0;
};