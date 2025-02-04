#include "sfuns.h"

double fgamma(double x){
    if(x<0) return std::numbers::pi/std::sin(std::numbers::pi*x)/fgamma(1-x); // Euler's reflection formula
    return std::exp(ln_gamma(x));
}

double ln_gamma(double x){
    if(x <= 0) return std::numbers::pi/std::sin(std::numbers::pi*x)/fgamma(1-x); // Euler's reflection formula
    if(x<9) return ln_gamma(x+1) - std::log(x); // Recurrence relation
    double value = x*std::log(x+1/(12*x-1/x/10))-x+std::log(2*std::numbers::pi/x)/2;
    return value;
}