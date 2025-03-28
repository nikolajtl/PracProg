#include <iostream>
#include <fstream>
#include <vector>
#include <cmath> // Contains NAN, nan and std::isnan
#include <functional> // For std::function
#include <limits> // For infinities
#include <utility> // For std::pair

double pi = 3.14159265358979323846;
double sqrt_pi = 1.7724538509055160272981674833411451827975494561223871282138077898529112845;
double pos_inf = std::numeric_limits<double>::infinity();
double neg_inf = -std::numeric_limits<double>::infinity();

std::vector<double> linspace(double start, double end, int num_points) {
    std::vector<double> result;
    double step = (end - start) / (num_points - 1);  // Calculate step size

    for (int i = 0; i < num_points; ++i) {
        result.push_back(start + i * step);
    }
    return result;
};



int Ncalls = 0;
std::pair<double, double> integrate (
    std::function<double(double)> f,
    double a, double b, double acc=1e-3, double eps=1e-3,
    double f2=NAN, double f3=NAN // NaN indicates first call
) {
    Ncalls++;

    double h=b-a;
    if(std::isnan(f2)){ f2=f(a+2*h/6); f3=f(a+4*h/6); } // first call, no points to reuse
    double f1 = f(a+h/6), f4 = f(a+5*h/6);
    double Q = (2*f1+f2+f3+2*f4)/6*(b-a); // higher order rule
    double q = (  f1+f2+f3+  f4)/4*(b-a); // lower order rule
    double err = std::abs(Q-q);
    if (err <= acc+eps*std::abs(Q)) return {Q, err};
    else {
        std::pair<double, double> int1 = integrate(f, a, (a+b)/2, acc/std::sqrt(2), eps, f1, f2); // int_i = {Q_i, err_i}
        std::pair<double, double> int2 = integrate(f, (a+b)/2, b, acc/std::sqrt(2), eps, f3, f4);
        return {int1.first + int2.first, std::sqrt(int1.second*int1.second + int2.second*int2.second)};
    };
};

std::pair<double, double> CCint (
    std::function<double(double)> f,
    double a, double b, double acc=1e-3, double eps=1e-3
) {
    std::function<double(double)> ftrans = [f, a, b](double theta)
    {return f((a+b)/2 + (b-a)/2*std::cos(theta))*std::sin(theta)*(b-a)/2;};
    return integrate(ftrans, 0, pi, acc, eps);
};

std::pair<double, double> infint(
    std::function<double(double)> f,
    double a, double b, double acc=1e-3, double eps=1e-3
) {
    if (std::isinf(a) && std::isinf(b)) { // Might cause problems that signs are ignored
        std::function<double(double)> ftrans = [f](double t)
        {return f(t/(1-t*t)) * (1+t*t)/((1-t*t)*(1-t*t));};
        return CCint(ftrans, -1, 1, acc, eps);
    };
    if (!std::isinf(a) && std::isinf(b)) {
        std::function<double(double)> ftrans = [f, a](double t)
        {return f(a + (1-t)/t)/t/t;};
        return CCint(ftrans, 0, 1, acc, eps);
    };
    if (std::isinf(a) && !std::isinf(b)) {
        std::function<double(double)> ftrans = [f, b](double t)
        {return f(b - (1-t)/t)/t/t;};
        return CCint(ftrans, 0, 1, acc, eps);
    };
    if (!std::isinf(a) && !std::isinf(b)) {
        return CCint(f, a, b, acc, eps);
    };
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

// The two test functions in B):
double CCtest1(double x) {
    return 1/std::sqrt(x);
};

double CCtest2(double x) {
    return std::log(x)/std::sqrt(x);
};

// Other test function(s) for B):
double gauss(double x) {
    return std::exp(-x*x);
};


// The error function, needs renaming to avoid issues with the cmath one
double my_erf(double z, double acc=1e-3, double eps=1e-3) {
    if (z<0) return -my_erf(-z);
    if (0<=z && z<=1) return 2/std::sqrt(pi) * integrate([](double x) {return std::exp(-x*x);}, 0, z, acc, eps).first;
    if (z>1) return 1 - 2/std::sqrt(pi) * integrate([z](double x) {return std::exp(-(z+(1-x)/x)*(z+(1-x)/x))/x/x;}, 0, 1, acc, eps).first;
};



int main() {
    double acc = 1e-4, eps = 1e-4; // May fail if required much better
    double testval1 = integrate(test1, 0, 1, acc, eps).first;
    double testval2 = integrate(test2, 0, 1, acc, eps).first;
    double testval3 = integrate(test3, 0, 1, acc, eps).first;
    double testval4 = integrate(test4, 0, 1, acc, eps).first;

    double CCtestval1 = CCint(CCtest1, 0, 1, acc, eps).first;
    double CCtestval2 = CCint(CCtest2, 0, 1, acc, eps).first;

    // std::cout << "Ncalls so far: " << Ncalls << std::endl;
    Ncalls = 0;
    std::pair<double, double> inftest = infint(gauss, neg_inf, pos_inf, acc, eps);
    double inftestval = inftest.first;
    double inftesterr = inftest.second;
    
    std::cout << "Value of the first A) test integral (should be 2/3): " << testval1 << std::endl;
    std::cout << "Value of the second A) test integral (should be 2): " << testval2 << std::endl;
    std::cout << "Value of the third A) test integral (should be pi): " << testval3 << std::endl;
    std::cout << "Value of the fourth A) test integral (should be -4): " << testval4 << std::endl;

    std::cout << "Value of the first B) test integral (should be 2): " << CCtestval1 << std::endl;
    std::cout << "Value of the second B) test integral (should be -4): " << CCtestval2 << std::endl;

    std::cout << "Value and error of the Gaussian integral, should be sqrt(pi), around 1.77245: "
    << inftestval << ", " << inftesterr << std::endl;
    std::cout << "Ncalls for Gaussian integral: " << Ncalls << std::endl; // Careful with the placing of the Ncalls reset
    std::cout << "Everything done for acc=" << acc << ", eps=" << eps << std::endl;

    std::ofstream outfile("out.txt");
    outfile << "Everything done for acc=" << acc << ", eps=" << eps << std::endl;
    outfile << "Value of the first A) test integral (should be 2/3): " << testval1 << std::endl;
    outfile << "Value of the second A) test integral (should be 2): " << testval2 << std::endl;
    outfile << "Value of the third A) test integral (should be pi): " << testval3 << std::endl;
    outfile << "Value of the fourth A) test integral (should be -4): " << testval4 << std::endl;
    outfile << "Value of the first B) test integral (should be 2): " << CCtestval1 << std::endl;
    outfile << "Value of the second B) test integral (should be -4): " << CCtestval2 << std::endl;
    outfile << "Value and error of the Gaussian integral, should be sqrt(pi), around 1.77245: "
    << inftestval << ", " << inftesterr << std::endl;
    outfile << "Actual error: " << std::abs(inftestval - sqrt_pi) << std::endl;
    outfile << "Ncalls for Gaussian integral: " << Ncalls << std::endl;
    outfile.close();

    std::vector<double> xs = linspace(-2, 2, 500);
    std::ofstream erffile("erfdata.txt");
    for (int i=0; i<xs.size(); i++) {
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