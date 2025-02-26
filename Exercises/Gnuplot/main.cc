#include <iostream>
#include <fstream>
#include "sfuns.h"

std::vector<double> tabulated_xs = {
    0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 3.0, 3.5
};

std::vector<double> tabulated_erfs = {
    0, 0.022564575, 0.045111106, 0.067621594, 0.090078126, 0.112462916, 0.222702589, 
    0.328626759, 0.428392355, 0.520499878, 0.603856091, 0.677801194, 0.742100965, 
    0.796908212, 0.842700793, 0.880205070, 0.910313978, 0.934007945, 0.952285120, 
    0.966105146, 0.976348383, 0.983790459, 0.989090502, 0.992790429, 0.995322265, 
    0.997020533, 0.998137154, 0.998856823, 0.999311486, 0.999593048, 0.999977910, 0.999999257
};


std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result;
    if (num <= 0) return result; // Return empty vector if num is invalid
    if (num == 1) { 
        result.push_back(start);
        return result;
    }

    double step = (end - start) / (num - 1);
    for (int i = 0; i<num; ++i) {
        result.push_back(start + i * step);
    }
    return result;
}


void writeDataToFile(const std::vector<double>& xvals, const std::vector<double>& yvals, const std::string& filename) {
    // Open the file in write mode
    std::ofstream outfile(filename);

    if (!outfile) {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }

    // Loop through both vectors and write pairs of values to the file
    for (size_t i = 0; i < xvals.size(); ++i) {
        outfile << xvals[i] << " " << yvals[i] << std::endl;  // Write x and y as space-separated values
    }

    // Close the file
    outfile.close();
}


int main() {
    std::vector<double> xvals_erf = linspace (-4, 4, 1000);
    std::vector<double> yvals_erf;
    for (double x : xvals_erf) yvals_erf.push_back(erf(x));

    std::vector<double> xvals_gamma = linspace(-4, 4, 1000);
    std::vector<double> yvals_gamma;
    std::vector<double> yvals_lngamma;
    for (double x : xvals_gamma) yvals_gamma.push_back(fgamma(x));
    for (double x : xvals_gamma) yvals_lngamma.push_back(ln_gamma(x));

    writeDataToFile(xvals_erf, yvals_erf, "data_erf.txt");
    writeDataToFile(tabulated_xs, tabulated_erfs, "tabulated_data.txt");
    writeDataToFile(xvals_gamma, yvals_gamma, "data_gamma.txt");
    writeDataToFile(xvals_gamma, yvals_lngamma, "data_lngamma.txt");
    return 0;
}