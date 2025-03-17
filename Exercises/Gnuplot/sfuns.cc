#include "sfuns.h"

double fgamma(double x){
    if(x<0) return std::numbers::pi/std::sin(std::numbers::pi*x)/fgamma(1-x); // Euler's reflection formula
    return std::exp(ln_gamma(x));
}

double ln_gamma(double x){
    // if(x < 0) {
    //     if (fgamma(x) >= 0) return std::log(fgamma(x));
    //     else return std::log(-fgamma(x)); // Tror jeg
    // }
    if(x<9) return ln_gamma(x+1) - std::log(x); // Recurrence relation
    double value = x*std::log(x+1/(12*x-1/x/10))-x+std::log(2*std::numbers::pi/x)/2;
    return value;
}

double erf(double x){
    if(x<0) return -erf(-x);
    std::vector<double> a={0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429};
    double t=1/(1+0.3275911*x);
    double sum=t*(a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*a[4])))); /* the right thing */
    return 1-sum*std::exp(-x*x);
}