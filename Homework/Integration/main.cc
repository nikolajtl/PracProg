#include <iostream>
#include <fstream>
#include <vector>
#include <cmath> // Contains NAN, nan and std::isnan
#include <functional> // For std::function

double pi = 3.14159265358979323846;

std::vector<double> linspace(double start, double end, int num_points) {
    std::vector<double> result;
    double step = (end - start) / (num_points - 1);  // Calculate step size

    for (int i = 0; i < num_points; ++i) {
        result.push_back(start + i * step);
    }
    return result;
};



double integrate (
    std::function<double(double)> f,
    double a, double b, double acc=1e-3, double eps=1e-3,
    double f2=NAN, double f3=NAN // NaN indicates first call
) {
    double h=b-a;
    if(std::isnan(f2)){ f2=f(a+2*h/6); f3=f(a+4*h/6); } // first call, no points to reuse
    double f1 = f(a+h/6), f4 = f(a+5*h/6);
    double Q = (2*f1+f2+f3+2*f4)/6*(b-a); // higher order rule
    double q = (  f1+f2+f3+  f4)/4*(b-a); // lower order rule
    double err = std::abs(Q-q);
    if (err <= acc+eps*std::abs(Q)) return Q;
    else return integrate(f, a, (a+b)/2, acc/std::sqrt(2), eps, f1, f2)
        + integrate(f, (a+b)/2, b, acc/std::sqrt(2), eps, f3, f4);
};



// The four test functions in A):
double test1(double x) {
    return std::sqrt(x);
};

double test2(double x) {
    return 1/std::sqrt(x);
};

double test3(double x) {
    return 4*std::sqrt(1-x*x);
};

double test4(double x) {
    return std::log(x)/std::sqrt(x);
};


// The error function, needs renaming to avoid issues with the cmath one
double my_erf(double z, double acc=1e-3, double eps=1e-3) {
    if (z<0) return -my_erf(-z);
    if (0<=z && z<=1) return 2/std::sqrt(pi) * integrate([](double x) {return std::exp(-x*x);}, 0, z, acc, eps);
    if (z>1) return 1 - 2/std::sqrt(pi) * integrate([z](double x) {return std::exp(-(z+(1-x)/x)*(z+(1-x)/x))/x/x;}, 0, 1, acc, eps);
};



int main() {
    double testval1 = integrate(test1, 0, 1, 1e-8, 1e-8);
    double testval2 = integrate(test2, 0, 1, 1e-8, 1e-8);
    double testval3 = integrate(test3, 0, 1, 1e-8, 1e-8);
    double testval4 = integrate(test4, 0, 1, 1e-8, 1e-8);
    
    std::cout << "Value of the first test integral (should be 2/3): " << testval1 << std::endl;
    std::cout << "Value of the first test integral (should be 2): " << testval2 << std::endl;
    std::cout << "Value of the first test integral (should be pi): " << testval3 << std::endl;
    std::cout << "Value of the first test integral (should be -4): " << testval4 << std::endl;

    std::ofstream outfile("out.txt");
    outfile << "Value of the first test integral (should be 2/3): " << testval1 << std::endl;
    outfile << "Value of the second test integral (should be 2): " << testval2 << std::endl;
    outfile << "Value of the third test integral (should be pi): " << testval3 << std::endl;
    outfile << "Value of the fourth test integral (should be -4): " << testval4 << std::endl;
    outfile.close();

    std::vector<double> xs = linspace(-2, 2, 500);
    std::ofstream erffile("erfdata.txt");
    for (int i; i<xs.size(); i++) {
        erffile << xs[i] << " " << my_erf(xs[i]) << std::endl;
    };
    erffile.close();

    std::vector<double> accs = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9};
    std::ofstream convfile("convdata.txt");
    for (double acc : accs) {
        convfile << acc << " " << std::abs(my_erf(1, acc, 0) - 0.84270079294971486934) << std::endl;
    };

    return 0;
};