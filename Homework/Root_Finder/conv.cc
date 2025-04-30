#include <iostream>
#include <fstream>
#include <vector>
#include <cmath> // Contains NAN, nan and std::isnan
#include <functional> // For std::function
#include <utility> // For std::pair
#include <random>
#include <optional> // Required for std::optional
#include "matrix.h"
#include "ode.h"



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


pp::matrix jacobian(
    std::function<pp::vector(pp::vector)> f,
    pp::vector x,
    std::optional<pp::vector> fx = std::nullopt,
    std::optional<pp::vector> dx = std::nullopt
) {
    if (!dx) {
       dx = pp::vector(x.size());
       std::transform(x.begin(), x.end(), dx->begin(), [](double xi) {
           return std::abs(xi) * std::pow(2, -26);
       });
    };
    if (!fx) fx = f(x);
    pp::matrix J = pp::matrix(x.size(), x.size());
    for (int j = 0; j < x.size(); j++) { // The dereferencing is related to the std::optional container
        (*dx)[j] = std::abs(x[j])*std::pow(2,-26); // Dereference dx here
        x[j] += (*dx)[j];  // Dereference dx here
        pp::vector df = f(x) - (*fx);  // Dereference fx here
        for (int i = 0; i < x.size(); i++) {
            J(i, j) = df[i] / (*dx)[j];  // Dereference dx here
        }
        x[j] -= (*dx)[j];  // Dereference dx here
    };
    return J;
};


pp::vector newton (
	std::function<pp::vector(pp::vector)> f,
    pp::vector start,
    double acc=1e-3,
    std::optional<pp::vector> δx = std::nullopt // Optional δx-vector for calculation of Jacobian
) {
    pp::vector x = start;
    pp::vector fx = f(x), z, fz;
    do{ /* Newton's iterations */
        if(fx.norm() < acc) break; /* job done */
        pp::matrix J = jacobian(f, x, fx, δx);
        QR QRJ(J);
        pp::vector Dx = QRJ.solve(-1*fx); /* Newton's step */
        double λ = 1;
        double λmin = 1.0/128.0; // Don't know what the size of this should be
        do{ /* linesearch */
            z=x+λ*Dx;
            fz = f(z);
            if (fz.norm() < (1-λ/2)*fx.norm()) break;
            if (λ < λmin) break;
            λ/=2;
        } while(true);
        x=z; fx=fz;
    } while(true);
    return x;
};


double Mfunc (
    pp::vector E, double rmax=8.0, double rmin=0.1, double rstep=0.1, double acc=1e-3, double eps=1e-3
) {
    auto schrodeq = [E](double r, pp::vector v) -> pp::vector {
        double f = v[0];
        double fprime = v[1];
        double fprimeprime = -2.0 * f * (E[0] + 1.0 / r);
        return pp::vector({fprime, fprimeprime});
    };

    std::pair<double, double> rinterval = {rmin, rmax}; // Interval for r-values
    pp::vector init({rmin-rmin*rmin, 1-2*rmin}); // Initial conditions, schrodeq(rmin)
    auto [rs, fs] = ode::driver(schrodeq, rinterval, init, rstep, acc, eps);

    return fs.back()[0]; // Return the last value of f(r), i.e. f(rmax)
};


int main() {
    std::ofstream rmaxfile("rmax.txt");

    std::vector<double> rmaxs = linspace(2.0, 7.0, 6);
    for (int j = 0; j < rmaxs.size(); ++j) {
        double rmax=rmaxs[j], rmin=0.1, rstep=0.1, acc=1e-3, eps=1e-3; // Take default parameters and wrap M so as to make it 'single-variable'
        auto wrappedMfunc = [=](pp::vector Evec) -> pp::vector {
            double E = Evec[0];
            return pp::vector({Mfunc(pp::vector({E}), rmax)});
        };

        double E = newton(wrappedMfunc, pp::vector({-1}), acc)[0];
        // std::cout << E << std::endl;

        auto schrodeq = [E](double r, pp::vector v) -> pp::vector {
            double f = v[0];
            double fprime = v[1];
            double fprimeprime = -2.0 * f * (E + 1.0 / r);
            return pp::vector({fprime, fprimeprime});
        };
    
        std::pair<double, double> rinterval = {rmin, rmax}; // Interval for r-values
        pp::vector init({rmin-rmin*rmin, 1-2*rmin}); // Initial conditions, schrodeq(rmin)
        auto [rs, fs] = ode::driver(schrodeq, rinterval, init, rstep, acc, eps);

        for (int i = 0; i < rs.size(); ++i) {
            rmaxfile << rs[i] << " " << fs[i][0] << std::endl;
        };
        if (j < rmaxs.size() - 1) {
            rmaxfile << "\n\n";
        };
    };

    rmaxfile.close();


    std::ofstream rminfile("rmin.txt");

    std::vector<double> rmins = linspace(0.1, 0.5, 5);
    for (int j = 0; j < rmins.size(); ++j) {
        double rmin=rmins[j], rmax=8.0, rstep=0.1, acc=1e-3, eps=1e-3; // Take default parameters and wrap M so as to make it 'single-variable'
        auto wrappedMfunc = [=](pp::vector Evec) -> pp::vector {
            double E = Evec[0];
            return pp::vector({Mfunc(pp::vector({E}), rmax, rmin)});
        };

        double E = newton(wrappedMfunc, pp::vector({-1}), acc)[0];
        // std::cout << E << std::endl;

        auto schrodeq = [E](double r, pp::vector v) -> pp::vector {
            double f = v[0];
            double fprime = v[1];
            double fprimeprime = -2.0 * f * (E + 1.0 / r);
            return pp::vector({fprime, fprimeprime});
        };
    
        std::pair<double, double> rinterval = {rmin, rmax}; // Interval for r-values
        pp::vector init({rmin-rmin*rmin, 1-2*rmin}); // Initial conditions, schrodeq(rmin)
        auto [rs, fs] = ode::driver(schrodeq, rinterval, init, rstep, acc, eps);

        for (int i = 0; i < rs.size(); ++i) {
            rminfile << rs[i] << " " << fs[i][0] << std::endl;
        };
        if (j < rmins.size() - 1) {
            rminfile << "\n\n";
        };
    };

    rminfile.close();


    std::ofstream accepsfile("acceps.txt");

    std::vector<double> acceps = {1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5};
    for (int j = 0; j < acceps.size(); ++j) {
        double rmin=0.1, rmax=8.0, rstep=0.1, acc=acceps[j], eps=acceps[j]; // Take default parameters and wrap M so as to make it 'single-variable'
        auto wrappedMfunc = [=](pp::vector Evec) -> pp::vector {
            double E = Evec[0];
            return pp::vector({Mfunc(pp::vector({E}), rmax, rmin, rstep, acc, eps)});
        };

        double E = newton(wrappedMfunc, pp::vector({-1}), acc)[0];
        // std::cout << E << std::endl;

        auto schrodeq = [E](double r, pp::vector v) -> pp::vector {
            double f = v[0];
            double fprime = v[1];
            double fprimeprime = -2.0 * f * (E + 1.0 / r);
            return pp::vector({fprime, fprimeprime});
        };
    
        std::pair<double, double> rinterval = {rmin, rmax}; // Interval for r-values
        pp::vector init({rmin-rmin*rmin, 1-2*rmin}); // Initial conditions, schrodeq(rmin)
        auto [rs, fs] = ode::driver(schrodeq, rinterval, init, rstep, acc, eps);

        for (int i = 0; i < rs.size(); ++i) {
            accepsfile << rs[i] << " " << fs[i][0] << std::endl;
        };
        if (j < acceps.size() - 1) {
            accepsfile << "\n\n";
        };
    };

    accepsfile.close();

    return 0;
};
