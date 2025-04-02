#include <iostream>
#include <fstream>
#include <vector>
#include <cmath> // Contains NAN, nan and std::isnan
#include <functional> // For std::function
#include <utility> // For std::pair
#include <random>
#include <cassert>
#include "matrix.h"

double pi = 3.1415926535897932384626433832795028841971693993751;

std::pair<double, double> plainmc(
    std::function<double(pp::vector)> f,
    pp::vector a, pp::vector b, int N
){
    int dim = a.size();
    double V=1;

    // Calculate the volume of the encompassing hypercube by multiplying axis lengths
    for(int i=0; i<dim; i++) V*=b[i]-a[i];

    double sum=0, sum2=0;
    pp::vector x(dim);

    // Somehow generates a number between zero and one
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (int i=0; i<N; i++) {
        for (int k=0; k<dim; k++) {
            x[k] = a[k] + dist(gen) * (b[k]-a[k]); // Random walk around the region
        };
        double fx=f(x); sum+=fx; sum2+=fx*fx;
    };

    double mean = sum/N, sigma=std::sqrt(sum2/N - mean*mean);
    std::pair<double, double> result(mean*V, sigma*V/std::sqrt(N));
    return result;
};


void quasirnd(pp::vector& x, bool reset=false) {
    static int dim = 0, n = 0;
    static std::vector<long double> alpha;

    if (reset) {
        dim = x.data.size(); // Infer dimension from x
        n = 0;
        alpha.resize(dim);

        for (int i=0; i<dim; i++) {
            double intpart; // Needed to call std::modf, not actually used for anything
            alpha[i] = std::modf(std::sqrt(pi + i), &intpart);
        };
    } else {
        assert(static_cast<int>(x.data.size()) == dim);
        n++;

        for (int i=0; i<dim; i++) {
            double intpart;
            x.data[i] = std::modf(n*alpha[i], &intpart);
        };
    };
};


std::pair<double, double> quasirndmc (
    std::function<double(pp::vector)> f,
    pp::vector a, pp::vector b, int N
) {
    int dim = a.size();
    double V=1;

    // Calculate the volume of the encompassing hypercube by multiplying axis lengths
    for (int i=0; i<dim; i++) V*=b[i]-a[i];

    // Generate the dim-dimensional grid of quasirandom points
    std::vector<pp::vector> grid;
    pp::vector x(N);
    quasirnd(x, true); // Generate the seed the first time around
    for (int i=0; i<dim; i++) {
        pp::vector x(N);
        quasirnd(x);
        grid.push_back(x);
    };

    double sum = 0;
    for (int i=0; i<N; i++) {
        pp::vector point(dim);
        for (int j=0; j<dim; j++) {
            // Scale each dimension's point to the range [a[j], b[j]]
            point[j] = a[j] + grid[j].data[i] * (b[j] - a[j]);
        };
        // Evaluate the function f at the point
        sum += f(point);
        std::cout << sum << std::endl;
    };

    // Estimate the integral: Average the function values and multiply by volume
    double integral = sum / N * V;
    
    return {integral, V};
};



double unitcirc (pp::vector v) {
    double x = v[0], y = v[1];
    if (x*x + y*y <= 1) {
        return 1;
    } else {
        return 0;
    };
};


// The super singular integration test function
double singintval = 1.3932039296856768591842462603255*pi*pi*pi;
double singfunc (pp::vector v) {
    double x = v[0], y = v[1], z = v[2];
    return 1 / (1-std::cos(x)*std::cos(y)*std::cos(z));
};


int main() {
    // Unit circle test
    pp::vector a = {-1, -1};
    pp::vector b = {1, 1};
    int N = 1e3;
    std::pair<double, double> result = quasirndmc(unitcirc, a, b, N);
    std::cout << "Actual error: " << std::abs(result.first - pi) << " , calculated error: " << result.second << std::endl;

    // Singular integral test
    // pp::vector a = {0, 0, 0};
    // pp::vector b = {pi, pi, pi};
    // int N = 1e3;
    // std::pair<double, double> result = quasirndmc(singfunc, a, b, N);
    // std::cout << "Value of singular integral (note the increase of a factor pi^3): " << result.first << ", expected value of integral: " << singintval << std::endl;
    // std::cout << "Actual error: " << std::abs(result.first - singintval) << " , calculated error: " << result.second << std::endl;

    // quasirnd(b, true);
    // b.print();
    // quasirnd(b);
    // b.print();
    // quasirnd(b);
    // b.print();

    return 0;
};