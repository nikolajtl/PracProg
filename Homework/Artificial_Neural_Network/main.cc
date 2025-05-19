#include <iostream>
#include <fstream>
#include <cmath>
#include "minimizer.h"


std::vector<double> linspace(double start, double end, int num_points) {
    std::vector<double> result;
    double step = (end - start) / (num_points - 1);  // Calculate step size

    for (int i = 0; i < num_points; ++i) {
        result.push_back(start + i * step);
    }
    return result;
};



struct ann {
    int n; // Number of hidden neurons (3)
    std::function<double(double)> f;
    pp::vector params;
    pp::vector w, a, b;

    ann (int n) : n(n), f( [](double x) { return x*std::exp(-x*x); } ) { // Gaussian wavelet activation function included in the constructor
        // Optionally initialize params here, might want different defaults in future
        std::vector<double> p(3*n, 1.0);
        params = pp::vector (p);
        // params.print();

        for (int i=0; i<n; i++) {
            w.push_back(params[i]);
            a.push_back(params[i+n]);
            b.push_back(params[i+2*n]);
        };
        // w.print(); a.print(); b.print();
    };

    double response(double x, const pp::vector& p) {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            double wi = p[i];
            double ai = p[i + n];
            double bi = p[i + 2 * n];
            sum += wi * f((x - ai) / bi);
        };
        return sum;
    };

    double cost(std::vector<double> xs, std::vector<double> ys, pp::vector p) {
        double cost = 0;
        for (size_t k = 0; k < xs.size(); ++k) {
            double diff = response(xs[k], p) - ys[k];
            cost += diff * diff;
        };
        return cost;
    };

    // Note that the accuracy here functions as the PARAMETER distance to the optimum, NOT the cost functions distance to zero
    void train (std::vector<double> xs, std::vector<double> ys, pp::vector guess, double acc=1e-3) {
        // Implement training procedure, using minimization algorithm due to laziness
        std::function<double(pp::vector)> costwrap = [&](pp::vector p) {
            return cost(xs, ys, p);
        };

        std::pair<pp::vector, int> min = min::newton(costwrap, guess, acc);
        params = min.first;
        std::cout << "Number of newton steps to train network for n=" << n << " neurons, acc=" << acc << "; " << min.second << std::endl;
        update_weights_from_params();
    };

    // Helper function that actually updates the parameters
    void update_weights_from_params() {
        w = pp::vector(); // Clear old weights
        a = pp::vector();
        b = pp::vector();
        for (int i = 0; i < n; ++i) {
            w.push_back(params[i]);
            a.push_back(params[i + n]);
            b.push_back(params[i + 2*n]);
        };
    };


    double response_int (double x, const pp::vector& p) {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            double wi = p[i];
            double ai = p[i + n];
            double bi = p[i + 2 * n];
            double ui = (x - ai) / bi;
            sum += wi * bi * std::exp(-ui*ui); // Analytical expression; must be re-done if activation function is changed
        };
        return -0.5*sum;
    };

    double response_deriv (double x, const pp::vector& p) {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            double wi = p[i];
            double ai = p[i + n];
            double bi = p[i + 2 * n];
            double ui = (x - ai) / bi;
            sum += (wi / bi) * (1 - 2*ui*ui) * std::exp(-ui*ui); // Analytical expression; must be re-done if activation function is changed
        };
        return sum;
    };

    double response_2deriv (double x, const pp::vector& p) {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            double wi = p[i];
            double ai = p[i + n];
            double bi = p[i + 2 * n];
            double ui = (x - ai) / bi;
            sum += (wi / (bi*bi)) * (4*ui*ui*ui - 6*ui) * std::exp(-ui*ui); // Analytical expression; must be re-done if activation function is changed
        };
        return sum;
    };
};


// The example function given in A)
std::function<double(double)> g = [](double x) { return std::cos(5*x - 1) * std::exp(-x*x); };


int main () {
    std::vector<double> xs = linspace(-1.5, 1.5, 200);
    std::vector<double> ys;
    for (double x : xs) {
        ys.push_back(g(x));
    };

    int n = 3; // Number of neurons
    ann neunet(n);

    pp::vector guess = {
        -1, 1.5, -0.7, // w
        -0.7, -0.2, 0.5, // a
        0.5, 0.5, 0.5  // b
    };

    std::cout << "Cost function value at starting guess: " << neunet.cost(xs, ys, guess) << std::endl;
    neunet.train(xs, ys, guess, 1);
    // neunet.train(xs, ys, neunet.params, 1); // Use default parameters as the guess (probably not great)
    std::cout << "Cost function value after training: " << neunet.cost(xs, ys, neunet.params) << std::endl;
    // neunet.params.print();

    std::ofstream valfile("values.txt");
    for (int i=0; i<xs.size(); i++) {
        valfile << xs[i] << " " << ys[i]
        << " " << neunet.response(xs[i], neunet.params)
        << " " << neunet.response_int(xs[i], neunet.params)
        << " " << neunet.response_deriv(xs[i], neunet.params)
        << " " << neunet.response_2deriv(xs[i], neunet.params)
        << std::endl;

        // Plot the guess for debugging:
        // valfile << xs[i] << " " << ys[i] << " " << neunet.response(xs[i], guess) << std::endl;
    };


    return 0;
};