#include <cmath>
#include <numbers>
#include <iostream>
#include "sfuns.h"

int main(){
    double sqrt2 = std::sqrt(2.0);
    double fifthroot2 = std::pow(2.0, 1.0/5.0);
    double e_to_pi = std::exp(std::numbers::pi);
    double pi_to_e = std::pow(std::numbers::pi, std::numbers::e);
    std::cout << "sqrt(2)^2 = " << sqrt2*sqrt2 << std::endl;
    std::cout << "(2^(1/5))^5 = " << std::pow(fifthroot2, 5.0) << std::endl;
    std::cout << "pi = " << std::log(e_to_pi) << std::endl;
    std::cout << "pi = " << std::pow(pi_to_e, 1/std::numbers::e) << std::endl;
    std::cout << "The actual calculated values above: " << sqrt2 << " " << fifthroot2 << " "<< e_to_pi << " "<< pi_to_e << " " << std::endl;
    std::cout << "Values of Gamma(n) from n=1 to n=10:" << std::endl;
    for (int n = 1; n <= 10; n++){
        std::cout << fgamma(n) << std::endl;
    }
}