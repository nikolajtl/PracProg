#include <iostream>
#include <complex>
#include <numbers>

bool approx(const std::complex<double>& a, const std::complex<double>& b, double abs_prec = 1e-9, double rel_prec = 1e-9) {
    // Compare the real parts
    if (std::abs(std::real(b) - std::real(a)) <= abs_prec) {
        // Compare the imaginary parts
        if (std::abs(std::imag(b) - std::imag(a)) <= abs_prec) return true;
        if (std::abs(std::imag(b) - std::imag(a)) / std::max(std::abs(std::imag(a)), std::abs(std::imag(b))) <= rel_prec) return true;
    };
    if (std::abs(std::real(b) - std::real(a)) / std::max(std::abs(std::real(a)), std::abs(std::real(b))) <= rel_prec) {
        // Compare the imaginary parts
        if (std::abs(std::imag(b) - std::imag(a)) <= abs_prec) return true;
        if (std::abs(std::imag(b) - std::imag(a)) / std::max(std::abs(std::imag(a)), std::abs(std::imag(b))) <= rel_prec) return true;
    }

    return false;
}

std::complex<double> plus_i(0, 1);
std::complex<double> minus_i(0, -1);

std::complex<double> c_pi(std::numbers::pi, 0);

std::complex<double> root_neg1 = sqrt( std::complex<double>(-1, 0) );
std::complex<double> root_i = sqrt(plus_i);
std::complex<double> exp_i = exp(plus_i);
std::complex<double> exp_i_pi = exp(plus_i*c_pi);
std::complex<double> i_to_i = pow(plus_i, plus_i);
std::complex<double> log_i = log(plus_i);
std::complex<double> sin_i_pi = sin(plus_i*c_pi);

int main() {
    std::cout << "Square root of -1: " << root_neg1
    << ", it seems to pick the 'positive' version: " << approx(root_neg1, minus_i) << std::endl;
    std::cout << "ln(i): " << log_i << std::endl;
    std::cout << "Square root of i: " << root_i << ", the reference states that the values are chosen in the right half-plane." <<std::endl;
    std::cout << "i^i: " << i_to_i << std::endl;
    return 0;
}