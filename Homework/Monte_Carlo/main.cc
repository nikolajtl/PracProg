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
double e = 2.71828182845904523536028747135266249775724709369995957;

std::vector<int> prime_numbers(int n) {
    std::vector<int> primes;
    int candidate = 2;
    while (primes.size() < n) {
        bool candidate_is_prime = true;
        for (int p : primes) {
            if (p * p > candidate) break;
            if (candidate % p == 0) {
                candidate_is_prime = false;
                break;
            };
        };
        if (candidate_is_prime) {
            primes.push_back(candidate);
        };
        candidate++;
    };
    return primes;
};


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


void quasirnd_irr(pp::vector& x, double irrnr, bool reset=false) {
    static int dim = 0, n = 0;
    static std::vector<long double> alpha;

    if (reset) {
        dim = x.data.size(); // Infer dimension from x
        n = 0;
        alpha.resize(dim);

        for (int i=0; i<dim; i++) {
            double intpart; // Needed to call std::modf, not actually used for anything
            alpha[i] = std::modf(std::sqrt(irrnr + i), &intpart);
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

void quasirnd_prime(pp::vector& x, int offset=0, bool reset=false) {
    static int dim = 0, n = 0;
    static std::vector<long double> alpha;

    if (reset) {
        dim = x.data.size(); // Infer dimension from x
        n = 0;
        alpha.resize(dim);
        std::vector<int> primes = prime_numbers(dim+offset);

        for (int i=0; i<dim; i++) {
            double intpart; // Needed to call std::modf, not actually used for anything
            alpha[i] = std::modf(std::sqrt(primes[i+offset]), &intpart);
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

    pp::vector q(dim);
    quasirnd_prime(q, 0, true); // Generate the seed the first time around
    pp::vector x(dim);
    
    double sum1 = 0;
    for (int i=0; i<N; i++) {
        for (int j=0; j<dim; j++) {
            quasirnd_prime(q);
            x[j] = a[j] + q[j] * (b[j]-a[j]);
        };
        // Evaluate the function f at the point
        sum1 += f(x);
    };

    // Estimate the integral: Average the function values and multiply by volume
    double integral1 = sum1 / N * V;

    quasirnd_prime(q, 5, true); // Generate the seed the first time around

    double sum2= 0;
    for (int i=0; i<N; i++) {
        for (int j=0; j<dim; j++) {
            quasirnd_prime(q, 5); // Offset only necesarry for reset
            x[j] = a[j] + q[j] * (b[j]-a[j]);
        };
        // Evaluate the function f at the point
        sum2 += f(x);
    };

    double integral2 = sum2 / N * V;
    double err = std::abs(integral1 - integral2);
    
    return {integral1, err};
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
    // pp::vector a = {-1, -1};
    // pp::vector b = {1, 1};
    // int N = 1e5;
    // std::pair<double, double> result = quasirndmc(unitcirc, a, b, N);
    // std::cout << "Actual error: " << std::abs(result.first - pi) << " , calculated error: " << result.second << std::endl;

    // Singular integral test
    pp::vector a = {0, 0, 0};
    pp::vector b = {pi, pi, pi};
    std::vector<int> Nvals = {10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 10000000, 20000000};
    std::ofstream datafile("data.txt");
    for (int N : Nvals) {
        std::pair<double, double> result_plain = plainmc(singfunc, a, b, N);
        std::pair<double, double> result_quasi = quasirndmc(singfunc, a, b, N);
        datafile << N << " " << std::abs(result_plain.first - singintval)/singintval
        << " " << result_plain.second/singintval << " " << result_quasi.second/singintval << std::endl;
    };
    datafile.close();

    int N = Nvals.back(); // Last entry
    std::pair<double, double> result_plain = plainmc(singfunc, a, b, N);
    std::cout << "Results for plain Monte Carlo:" << std::endl;
    std::cout << "Value of singular integral (note the increase of a factor pi^3): " << result_plain.first << ", expected value of integral: " << singintval << std::endl;
    std::cout << "Actual error: " << std::abs(result_plain.first - singintval) << " , calculated error: " << result_plain.second << std::endl;
    std::cout << std::endl;
    std::pair<double, double> result_quasi = quasirndmc(singfunc, a, b, N);
    std::cout << "Results for quasi-random Monte Carlo:" << std::endl;
    std::cout << "Value of singular integral (note the increase of a factor pi^3): " << result_quasi.first << ", expected value of integral: " << singintval << std::endl;
    std::cout << "Actual error: " << std::abs(result_quasi.first - singintval) << " , calculated error: " << result_quasi.second << std::endl;

    return 0;
};