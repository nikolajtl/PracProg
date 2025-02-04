#include <iostream>
#include <cmath>
#include <algorithm> // For std::max

int n=1, m=-2;
float eps_float = 1.0;
double eps_doub = 1.0;

bool approx(double a, double b, double abs_prec=1e-9, double rel_prec=1e-9){
        if (std::abs(b-a) <= abs_prec) return true;
        if (std::abs(b-a)/std::max(a, b) <= rel_prec) return true;
        return false;
}

int main(){
    while(n+1 > n){
        n += 1;
    }

    while(m-1 < m){
        m -= 1;
    }

    while(1.0 + eps_float != 1.0){
        eps_float /= 2;
    }

    while(1.0 + eps_doub != 1.0){
        eps_doub /= 2;
    }

    std::cout << "Maximum representable integer: " << n << std::endl;
    std::cout << "Minimum representable integer: " <<  m << std::endl;
    std::cout << "Float epsilon: " << eps_float << std::endl;
    std::cout << "Double epsilon: " << eps_doub << std::endl;
    std::cout << "2^(-23): " << std::pow(2, -23) << std::endl;
    std::cout << "2^(-53): " << std::pow(2, -53) << std::endl;

    double tiny = eps_doub / 2;
    double a = 1 + tiny + tiny;
    double b = tiny + tiny + 1;

    std::cout << "a: " << a << ", b: " << b << std::endl;
    std::cout << "a == b: " << (a == b) << std::endl;
    std::cout << "a > 1: " << (a > 1) << std::endl;
    std::cout << "b > 1: " << (b > 1) << std::endl;

    double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
    double d2 = 8*0.1;

    std::cout << "d1, d2 'strictly' equal? " << (d1 == d2) << std::endl;
    std::cout << "d1, d2 approximately equal? " << approx(d1, d2) << std::endl;
}