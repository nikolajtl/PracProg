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

    ann (int n) : n(n), f( [](double x) { return std::exp(-x*x); } ) { // Gaussian activation function included in the constructor
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

    // Note that the phi integral is included in its function definition, so it's just called here
    double cost(std::vector<double> xs, double c, double y_c, double y_c_prime,
                std::function<double(std::function<double(double)>, std::vector<double>)> phi, pp::vector p,
                double delta_c = 0.01 // The x-step used for calculating the derivative at x=c
            ) {
        double cost = 0;
        double alpha = 100, beta = 100; // Initial condition cost weight parameters
        // for (size_t k = 0; k < xs.size(); ++k) {
        //     double phival = phi(response(xs[k], p));
        //     cost += phival;
        // };
        std::function<double(double)> response_wrap = [this, p](double x) {
            return this->response(x, p);
        };
        cost += phi(response_wrap, xs);
        cost += alpha * (response(c, p) - y_c)*(response(c, p) - y_c);
        double Fp_prime_c = (response(c+delta_c, p) - response(c-delta_c, p)) / (2*delta_c);
        cost += beta * (Fp_prime_c - y_c_prime)*(Fp_prime_c - y_c_prime);
        return cost;
    };

    // Note that the accuracy here functions as the PARAMETER distance to the optimum, NOT the cost functions distance to zero
    void train (std::vector<double> xs, double c, double y_c, double y_c_prime,
                std::function<double(std::function<double(double)>, std::vector<double>)> phi,
                pp::vector guess, double acc=1e-3
            ) {
        // Implement training procedure, using minimization algorithm due to laziness
        std::function<double(pp::vector)> costwrap = [&](pp::vector p) {
            return cost(xs, c, y_c, y_c_prime, phi, p);
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

};


double omega = 5.0, eta = 1.3;
// Dampened harmonic oscillator differential equation, eta=0 for undamped
double damp_harm_diffeq(std::function<double(double)> y, std::vector<double> xs) {
    double phi_sq_integral = 0.0;

    std::vector<double> ys;
    std::vector<double> y_primes;
    std::vector<double> y_2primes;

    for (int i = 1; i < xs.size() - 1; ++i) {
        double h1 = xs[i] - xs[i - 1];
        double h2 = xs[i + 1] - xs[i];
        
        double y_prev = y(xs[i - 1]);
        double y_curr = y(xs[i]);
        double y_next = y(xs[i + 1]);
        
        double y_prime = (y_next - y_prev) / (h1 + h2);
        double y_2prime = (y_next - 2 * y_curr + y_prev) / ((0.5 * (h1 + h2)) * (0.5 * (h1 + h2)));
        
        double phi = y_2prime + eta * y_prime + omega * omega * y_curr;
        double dx = (xs[i + 1] - xs[i - 1]) / 2.0;
        phi_sq_integral += phi * phi * dx;
    };

    // for (int i=0; i<xs.size()-2; i++) { // Remove the last two to avoid out-of-bounds issues
    //     double xi = xs[i], xip1 = xs[i+1], xip2 = xs[i+2];
    //     double yi = y(xi), yip1 = y(xip1), yip2 = y(xip2);
        
    //     double y_prime_i = (yip1 - yi) / (xip1 - xi);
    //     double y_prime_ip1 = (yip2 - yip1) / (xip2 - xip1);

    //     double y_2prime_i = (y_prime_ip1 - y_prime_i) / (xip1 - xi);

    //     ys.push_back(yi);
    //     y_primes.push_back(y_prime_i);
    //     y_2primes.push_back(y_2prime_i);

    //     double phi = y_2prime_i + eta * y_prime_i + omega*omega*yi;
    //     double phi_sq = phi*phi;
    //     phi_sq_integral += phi_sq * (xip1 - xi);

    //     // phi += (y_2prime_i + eta * y_prime_i + omega*omega*yi) * xi; // Estimate the phi by 'integrating' over x
    // };

    return phi_sq_integral;
};


double analytic_sol (double x) {
    return std::exp(-0.5*eta*x)*std::cos(std::sqrt(omega*omega - (eta*eta)/4)*x);
};


double sign_int (int n) {return (n % 2 == 0) ? 1.0 : -1.0;};


int main() {
    std::vector<double> xs = linspace(0, 3, 100);
    double c=0, y_c=1, y_c_prime=0;

    int n = 8; // Number of neurons
    ann neunet(n);

    pp::vector guess; // = {
        // 1, -0.95, 0.9, -0.85, 0.8, -0.75, 0.7, -0.65, 0.6, -0.55, 0.5, -0.45, 0.4, -0.35, 0.3, -0.25, 0.2, -0.15, 0.1, 0.05 // w
        // -0.7, -0.2, 0.5, // a
        // 0.5, 0.5, 0.5  // b
    // };
    for (int i = 0; i < n; ++i) guess.push_back(sign_int(i)*0.8*(1 - 1.2*(double)i/n)); // Weights w_i
    // for (int i = 0; i < n; ++i) guess.push_back(0.1*sign_int(i));           // weights w_i
    // for (int i = 0; i < n; ++i) guess.push_back(0.2);
    for (int i = 0; i < n; ++i) guess.push_back(3.0*i/(n-1)); // centers a_i spaced in [0.0, 3.0]
    for (int i = 0; i < n; ++i) guess.push_back(0.2);           // widths b_i

    // pp::vector guess = {
    //     // Weights (w)
    //     -0.5082370343271161, -0.7263892530794906, 0.35451263317917436, -1.0873351246905898,
    //     0.6838019767908963, -1.076440426574263, 0.6308405345517706, -0.7859158936039656,
    //     0.33519658771813454, -0.36993535249082343, -0.03048824025848771, 0.013733253032212822,
    //     -0.32249331301648243, 0.25624405231595515, -0.4634658481911347, 0.3199290092258974,
    //     -0.4488264018173084, 0.2309462487799942, -0.3251169349460537, 0.022608949331951494,
        
    //     // Centers (a)
    //     0.0, 0.15789473684210525, 0.3157894736842105, 0.47368421052631576, 0.631578947368421,
    //     0.7894736842105263, 0.9473684210526315, 1.1052631578947367, 1.263157894736842,
    //     1.4210526315789473, 1.5789473684210527, 1.7368421052631577, 1.894736842105263,
    //     2.052631578947368, 2.2105263157894735, 2.3684210526315788, 2.526315789473684,
    //     2.6842105263157894, 2.8421052631578947, 3.0,
        
    //     // Widths (b)
    //     0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
    //     0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
    // };

    // for (int i = 0; i < n; ++i) guess.push_back(0.08);           // widths b_i

    std::cout << "Cost function value at starting guess: "
              << neunet.cost(xs, c, y_c, y_c_prime, damp_harm_diffeq, guess) << std::endl;
    // neunet.train(xs, c, y_c, y_c_prime, damp_harm_diffeq, neunet.params, 1e4); // Use default parameters as the guess (probably not great)
    neunet.train(xs, c, y_c, y_c_prime, damp_harm_diffeq, guess, 5);
    std::cout << "Cost function value after training: "
              << neunet.cost(xs, c, y_c, y_c_prime, damp_harm_diffeq, neunet.params) << std::endl;
    
    std::ofstream valfile("diffeq_values.txt");
    for (int i=0; i<xs.size(); i++) {
        valfile << xs[i] << " " << neunet.response(xs[i], neunet.params) << " " << analytic_sol(xs[i]) << std::endl;

        // Plot the guess for debugging:
        // valfile << xs[i] << " " << neunet.response(xs[i], guess) << " " << analytic_sol(xs[i]) << std::endl;
    };

    return 0;
};