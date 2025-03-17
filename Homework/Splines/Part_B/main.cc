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


double quadterp(double z, std::vector<double> xs, std::vector<double> ys) { // Produces the value val=f_linspline(z) given the lists of (x_i, y_i) points
    int i=binsearch(z, xs);
    double dx = xs[i+1] - xs[i];
    double dy = ys[i+1] - ys[i];
    double p = dy/dx;
    double c=0;
    for (int j=0; j<i; j++) {
        double dx_j = dx;
        double dx_jp1 = xs[j+2] - xs[j+1]; // jp1 read as "j plus one"
        double p_j = p;
        double dy_jp1 = (ys[j+2] - ys[j+1]) / dx_jp1;
        double p_jp1 = dy_jp1 / dx_jp1;
        c = (p_jp1 - p_j - c*dx_j) / dx_jp1;
    };
    return ys[i]+p*(z-xs[i]) + c*(z-xs[i])*(z-xs[i+1]);
};


double quadtegral(double z, std::vector<double> xs, std::vector<double> ys) {
    int i=binsearch(z, xs);
    if (z==xs[0]) return 0; // Bad debugging
    double integral = 0;
    for (int j=0; j<i; j++) {
        double dx_j = xs[j+1] - xs[j];
        double dy_j = ys[j+1] - ys[j];
        double p_j = dy_j/dx_j;
        // Calculate c_j:
        double c_j = 0;
        for (int k=0; k<j; k++) {
            double dx_k = xs[k+1] - xs[k];
            double dx_kp1 = xs[k+2] - xs[k+1];
            double dy_k = ys[k+1] - ys[k];
            double dy_kp1 = ys[k+2] - ys[k+1];
            double p_k = dy_k / dx_k;
            double p_kp1 = dy_kp1 / dx_kp1;
            c_j = (p_kp1 - p_k - c_j*dx_k) / dx_kp1;
        };
        double b_j = p_j - c_j*dx_j;
        integral += ys[j]*dx_j + b_j*dx_j*dx_j / 2 + c_j*dx_j*dx_j*dx_j / 3; // Eqn. 1.15
    };
    double dx = z - xs[i];
    double dy = quadterp(z, xs, ys) - ys[i];
    double p = dy/dx;
    double c=0;
    for (int k=0; k<i; k++) {
        double dx_k = xs[k+1] - xs[k];
        double dx_kp1 = xs[k+2] - xs[k+1];
        double dy_k = ys[k+1] - ys[k];
        double dy_kp1 = ys[k+2] - ys[k+1];
        double p_k = dy_k / dx_k;
        double p_kp1 = dy_kp1 / dx_kp1;
        c = (p_kp1 - p_k - c*dx_k) / dx_kp1;
    };
    double b = p - c*dx;
    integral += ys[i]*dx + b*dx*dx / 2 + c*dx*dx*dx / 3;
    return integral;
};


double spline(double x, std::vector<double> xs, std::vector<double> ys) {
    double y = quadterp(x, xs, ys);
    return y;
};


int main() {
    // Dummy x-vals
    std::vector<double> xs = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    // Make sin(x) points
    std::vector<double> ys(xs.size());
    std::transform(xs.begin(), xs.end(), ys.begin(), [](double n) { return std::sin(n); });


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
                   [&](double x) { return quadtegral(x, xs, ys); });
    
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

    std::cout << quadterp(0, xs, ys) << quadtegral(0, xs, ys) << std::endl;

    return 0;
};