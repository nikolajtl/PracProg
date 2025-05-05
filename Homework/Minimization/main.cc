#include <iostream>
#include <fstream>
#include <vector>
#include <cmath> // Contains NAN, nan and std::isnan
#include <functional> // For std::function
#include <utility>
#include "matrix.h"


struct QR {
    pp::matrix A, Q, R;

    QR(const pp::matrix& A) : A(A), Q(A.size1(), A.size2()), R(A.size2(), A.size2()) {
        // this->A.print();
        for (int i = 0; i < A.size2(); ++i) {
            pp::vector a = A.cols[i];

            for (int j = 0; j < i; ++j) {
                double r = a.dot(Q.cols[j]);
                R.set(j, i, r);
                a = a - (Q.cols[j] * r);
            }

            double norm = a.norm();
            if (norm == 0) throw std::runtime_error("Linear dependency detected");

            Q.cols[i] = a / norm;
            R.set(i, i, norm);
        }
    };

    pp::vector solve(const pp::vector& b) {
        pp::vector c = Q.transpose() * b;

        std::vector<double> x_vals(R.size1(), 0);
        for (int i = R.size1() - 1; i >= 0; --i) {
            double sum = 0;
            for (int j = i + 1; j < R.size1(); ++j) {
                sum += R.get(i, j) * x_vals[j];
            }
            x_vals[i] = (c.data[i] - sum) / R.get(i, i);

        }
        pp::vector sol(R.size1());
        sol.data = x_vals;
        return sol;
    };

    double det() {
        int n = R.size1();
        double det = 1;
        for (int i = 0; i < n; i++) {
            det *= R.get(i, i);
        }
        return det;
    };

    pp::matrix inverse() {
        if (det() == 0) throw std::runtime_error("Matrix is singular");
        if (A.size1() != A.size2()) throw std::runtime_error("Matrix is not square");

        pp::matrix inv(A.size2(), A.size1());
        for (int i = 0; i < A.size2(); i++) {
            // Generate the i-th unit vector
            pp::vector e_i;
            for (int j = 0; j < A.size1(); j++) {
                e_i.data.push_back(i == j ? 1 : 0);
            }

            // Solve A * x_i = e_i
            pp::vector x = solve(e_i);
            for (int j = 0; j < A.size1(); j++) {
                inv.set(j, i, x[j]); // x_i is the i'th column of the inverse
            }
        }
        return inv;
    };

    pp::matrix Rinverse() {
        if (R.size1() != R.size2()) throw std::runtime_error("R is not square");
    
        int n = R.size1();
        pp::matrix Rinv(n, n);
    
        // Iterate through columns to compute the inverse
        for (int i = n - 1; i >= 0; --i) {
            // Solve R * x = e_i, where e_i is the i-th unit vector
            std::vector<double> x_vals(n, 0);
            x_vals[i] = 1.0 / R.get(i, i);
    
            for (int j = i - 1; j >= 0; --j) {
                double sum = 0;
                for (int k = j + 1; k < n; ++k) {
                    sum += R.get(j, k) * x_vals[k];
                }
                x_vals[j] = -sum / R.get(j, j);
            }
    
            for (int j = 0; j < n; ++j) {
                Rinv.set(j, i, x_vals[j]);
            }
        }
    
        return Rinv;
    };    
};

std::vector<double> linspace(double start, double end, int num_points) {
    std::vector<double> result;
    double step = (end - start) / (num_points - 1);  // Calculate step size

    for (int i = 0; i < num_points; ++i) {
        result.push_back(start + i * step);
    }
    return result;
};


// Gradient backup
// pp::vector gradient(std::function<double(pp::vector)> φ, pp::vector x) {
//     double φx = φ(x);
//     pp::vector gφ = pp::vector(x.size()); // Gradient of phi (construct as empty)
//     for (int i=0; i < x.size(); i++) {
//         double dxi = std::abs(x[i]) * std::pow(2, -26);
//         x[i] += dxi;
//         gφ[i] = (φ(x) - φx) / dxi;
//         x[i] -= dxi;
//     };
//     return gφ;
// };

// Hessian backup
// pp::matrix hessian(std::function<double(pp::vector)> φ, pp::vector x) {
//     pp::matrix H(x.size(), x.size());
//     pp::vector gφx = gradient(φ, x);
//     for (int j=0; j < x.size(); j++) {
//         double dxj = std::abs(x[j]) * std::pow(2, -17);
//         x[j] += dxj;
//         pp::vector dgφ = gradient(φ, x) - gφx;
//         for (int i=0; i < x.size(); i++) H(i, j) = dgφ[i] / dxj;
//         x[j] -= dxj;
//     };
//     return H;
// };


pp::vector gradient(std::function<double(pp::vector)> φ, pp::vector x, int mode=0) {
    pp::vector x_start = x;
    double φx = φ(x);
    double φx_plus = φx, φx_minus = φx;
    pp::vector gφ_forward = pp::vector(x.size()); // Gradient of phi (construct as empty)
    pp::vector gφ_central = pp::vector(x.size());
    for (int i=0; i < x.size(); i++) {
        double dxi = std::abs(x[i]) * std::pow(2, -26);
        x[i] += dxi;
        φx_plus = φ(x);
        gφ_forward[i] = (φx_plus - φx) / dxi;
        x[i] -= 2*dxi;
        φx_minus = φ(x);
        gφ_central[i] = (φx_plus - φx_minus) / (2*dxi);
        x[i] += dxi; // Reset for the next iteration of the loop
    };
    if (mode == 0) {
        return gφ_forward;
    } else if (mode == 1) {
        return gφ_central;
    } else {
        throw std::runtime_error("Invalid mode in gradient call; must be 0 for forward or 1 for central.");
    };
};


pp::matrix hessian(std::function<double(pp::vector)> φ, pp::vector x, int mode=0) {
    pp::matrix H_forward(x.size(), x.size());
    pp::matrix H_central(x.size(), x.size());
    pp::vector gφx_forward = gradient(φ, x, 0), gφx_central = gradient(φ, x, 1);
    pp::vector gφx_forward_plus = gφx_forward;
    pp::vector gφx_central_plus = gφx_central, gφx_central_minus = gφx_central;
    for (int j=0; j < x.size(); j++) {
        double dxj = std::abs(x[j]) * std::pow(2, -26);
        x[j] += dxj;
        gφx_forward_plus = gradient(φ, x, 0);
        gφx_central_plus = gradient(φ, x, 1);
        x[j] -= 2*dxj;
        gφx_central_minus = gradient(φ, x, 1);
        x[j] += dxj;

        pp::vector dgφ_forward = gφx_forward_plus - gφx_forward;
        pp::vector dgφ_central = gφx_central_plus - gφx_central_minus;

        for (int i=0; i < x.size(); i++) H_forward(i, j) = dgφ_forward[i] / dxj;
        for (int i=0; i < x.size(); i++) H_central(i, j) = dgφ_central[i] / (2*dxj);
    };
    if (mode == 0) {
        return H_forward;
    } else if (mode == 1) {
        return H_central;
    } else {
        throw std::runtime_error("Invalid mode in hessian call; must be 0 for forward or 1 for central.");
    };
};


std::pair<pp::vector, int> newton(std::function<double(pp::vector)> φ, pp::vector x, double acc=1e-3, int mode=0) {
    int steps = 0;
    while (true) {
        pp::vector g = gradient(φ, x, mode);

        if (g.norm() < acc) break; // Finished condition

        pp::matrix H = hessian(φ, x, mode);
        QR qrH(H);
        pp::vector dx = qrH.solve(-1*g); // dx not necesarilly small before mult. by λ

        double λ = 1.0;
        while (λ > 1.0/128.0) { // Backtrack linesearch
            if ( φ(x + λ*dx) < φ(x) ) break; // 'Good' step
            λ /= 2;
        };
        x = x + λ*dx;
        steps++;
    };
    return {x, steps};
};

void load_higgs_data (const std::string& filename,
                std::vector<double>& E,
                std::vector<double>& sigma,
                std::vector<double>& delta_sigma
                ) {
    std::ifstream infile(filename);
    if (!infile) throw std::runtime_error("Could not open file: " + filename);
    
    double e, s, d;
    while (infile >> e >> s >> d) {
        E.push_back(e);
        sigma.push_back(s);
        delta_sigma.push_back(d);
    };
};


// Rosenbrock's valley function
double rbvalley (pp::vector v) {
    double x = v[0];
    double y = v[1];
    double val = (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
    return val;
};

// Himmelblau's function
double himblau (pp::vector v) {
    double x = v[0];
    double y = v[1];
    double val = (x*x+y-11)*(x*x+y-11) + (x+y*y-7)*(x+y*y-7);
    return val;
};

// Breit-Wigner with parameters
double BW (double E, pp::vector params) {
    double m = params[0], Γ = params[1], A = params[2]; // Params order: m, Γ, A
    double val = A / ( (E-m)*(E-m) + Γ*Γ/4 );
    return val;
};

// χ^2 w.r.t. the data set
double χ2 (std::vector<double> E,
            std::vector<double> sigma,
            std::vector<double> sigma_err,
            pp::vector params) {
    // double m = params[0], Γ = params[1], A = params[2]; // Params order: m, Γ, A
    double chisq = 0;
    for (int i=0; i<E.size(); i++) {
        chisq += ( (sigma[i] - BW(E[i], params)) / sigma_err[i] ) * ( (sigma[i] - BW(E[i], params)) / sigma_err[i] );
    };
    return chisq;
};



int main() {
    double acc = 1e-6;
    std::vector<int> modes = {0, 1};
    
    pp::vector higgs_params;
    for (int mode : modes) {
        if (mode == 0) {
            std::cout << "Find all minima with forward difference: " << std::endl;
        };
        if (mode == 1) {
            std::cout << "Find all minima with central difference: " << std::endl;
        };

        std::pair<pp::vector, int> rbsol = newton(rbvalley, {4, 5}, acc, mode);
        std::pair<pp::vector, int> hbsol = newton(himblau, {2.5, 3.5}, acc, mode);
        pp::vector rbmin = rbsol.first, hbmin = hbsol.first;
        int rbsteps = rbsol.second, hbsteps = hbsol.second;

        std::cout << "Rosenbrock valley function minimum: " << std::endl;
        rbmin.print();
        std::cout << "Found in " << rbsteps << " steps from starting guess (4, 5)" << std::endl;

        std::cout << std::endl;

        std::cout << "One of the Himmelblau function's local minima: " << std::endl;
        hbmin.print();
        std::cout << "Found in " << hbsteps << " steps from starting guess (2.5, 3.5)" << std::endl;

        std::cout << std::endl;


        std::vector<double> E, sigma, sigma_err;
        load_higgs_data("higgs_data.txt", E, sigma, sigma_err);
        
        // Test if the data is loaded correctly (for debugging)
        // std::cout << "Loaded " << E.size() << " data points." << std::endl;
        // for (size_t i = 0; i < E.size(); ++i) {
        //     std::cout << "E = " << E[i]
        //               << ", σ = " << sigma[i]
        //               << ", Δσ = " << sigma_err[i] << std::endl;
        // };

        // Wrap chi2 to feed the minimization routine
        auto χ2_lambda = [&](pp::vector params) {
            return χ2(E, sigma, sigma_err, params);
        };

        std::pair<pp::vector, int> higgssol = newton(χ2_lambda, {125, 5, 25}, acc, mode);
        higgs_params = higgssol.first;
        int higgssteps = higgssol.second;

        std::cout << "Optimal parameters for the Higgs Breit-Wigner (m, Γ, A): " << std::endl;
        higgs_params.print();
        std::cout << "Found in " << higgssteps << " steps from starting guess (125, 5, 25). The χ2 value is "
                    << χ2(E, sigma, sigma_err, higgs_params) << std::endl;

        std::cout << std::endl;
    };
    std::cout << "All done for accuracy " << acc << ". It seems to vary from function to function which method is more effective." << std::endl;


    std::vector<double> Espace = linspace(100, 160, 500);
    std::ofstream fitfile("higgs_fit.txt");
    for (int i=0; i<Espace.size(); i++) {
        fitfile << Espace[i] << " " << BW(Espace[i], higgs_params) << std::endl;
        // fitfile << Espace[i] << " " << BW(Espace[i], {125, 5, 25}) << std::endl; // 'Fit by eye' for finding a good guess
    };
    fitfile.close();

    return 0;
};
