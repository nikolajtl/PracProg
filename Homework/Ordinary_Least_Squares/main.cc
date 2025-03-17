#include <iostream>
#include <fstream> // For writing out the data
#include <vector>
#include <cmath>
#include <random>
#include <functional> // For std::function
#include <algorithm> // For std::transform
#include "matrix.h"


struct QR {
    pp::matrix A, Q, R;

    QR(const pp::matrix& A) : A(A), Q(A.size1(), A.size2()), R(A.size2(), A.size2()) {
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


struct fit {
    pp::vector c; // Coefficients

    pp::matrix Sigma;
    pp::vector cerr;

    fit(const std::vector<std::function<double(double)>>& fs,
        std::vector<double> xs, std::vector<double> ys, std::vector<double> dys) {
        pp::matrix A; // Apparently not A() for empty constructor
        pp::vector b;

        // Construct A
        if (xs.size()==ys.size() && ys.size()==dys.size()) {
            for (int k=0; k<fs.size(); k++) {
                pp::vector A_k;
                for (int i=0; i<xs.size(); i++) {
                    double A_ki = fs[k](xs[i]) / dys[i];
                    A_k.push_back(A_ki);
                };
                A.push_back(A_k);
            };
        } else {
            throw std::runtime_error("Mismatch between lengths of the three data sets; x, y, dy");
        };

        // Construct b, and we've already checked that the sizes match
        for (int i=0; i<ys.size(); i++) {
            double b_i = ys[i] / dys[i];
            b.push_back(b_i);
        };

        QR QR_A(A);
        c = QR_A.solve(b);

        // Sigma = A^(-1) * A^(-1)^T = R^(-1) * R^(-1)^T
        pp::matrix Rinv = QR_A.Rinverse();
        pp::matrix RinvT = Rinv.transpose();
        Sigma = Rinv * RinvT;
        // Construct the vector of errors in the fitted values
        for (int i; i<Sigma.size1(); i++) {
            double err = std::sqrt( Sigma(i, i) );
            cerr.push_back(err);
        };
    };
};


double decay(double t, double ln_a, double lambda) {
    double a = std::exp(ln_a);
    return a * std::exp(-lambda * t);
};



int main() {
    std::vector<std::function<double(double)>> ln_y = {
        [](double t) { return 1.0; }, // ln_a
        [](double t) { return - t; }, // - lambda * t
        };
    
    std::vector<double> ts = {1, 2, 3, 4, 6, 9, 10, 13, 15};
    std::vector<double> Ns = {117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1};
    std::vector<double> dNs = {6, 5, 4, 4, 4, 3, 3, 2, 2};

    std::vector<double> ln_Ns(Ns.size());
    std::transform(Ns.begin(), Ns.end(), ln_Ns.begin(), [](double n) { return std::log(n); });

    std::vector<double> ln_dNs(dNs.size());
    std::transform(dNs.begin(), dNs.end(), ln_dNs.begin(), [](double n) { return std::log(n) / n; });

    fit expfit(ln_y, ts, ln_Ns, ln_dNs);
    std::cout << "Value of ln(a): " << expfit.c[0]
    << ". Error in ln(a): " << expfit.cerr[0]
    << ". Value of lambda: " << expfit.c[1]
    << ". Error in lambda: " << expfit.cerr[1]
    << ". We can turn this into a half-life with uncertainty by differentiating the expression T_(1/2)=ln(2)/lambda, wrt. lambda, and taking the absolute value (equiv. to squaring and taking square root). This yields T_(1/2) being "
    << std::log(2)/expfit.c[1]
    << " days, and the error being " << std::log(2)/(expfit.c[1]*expfit.c[1])*expfit.cerr[1] << " days. The modern value is around 3.6 days, so it fits nicely within the error."
    << std::endl;
    
    // Put it into a data file
    std::ofstream outFile("fit.txt");
    if (!outFile) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return 1;
    };

    double t_start = 0.0;
    double t_end = 16.0;
    double step = 0.1;  // Step size

    for (double t = t_start; t <= t_end; t += step) {
        double value = decay(t, expfit.c[0], expfit.c[1]);
        outFile << t << "\t" << value << "\n";
    };

    outFile.close();
    std::cout << "Wrote data for the fitted function to fit.txt" << std::endl;


    std::ofstream outFilepp("fit_ap_lp.txt");
    if (!outFilepp) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return 1;
    };

    for (double t = t_start; t <= t_end; t += step) {
        double value = decay(t, expfit.c[0]+expfit.cerr[0], expfit.c[1]+expfit.cerr[1]);
        outFilepp << t << "\t" << value << "\n";
    };

    outFilepp.close();
    std::cout << "Wrote data for the fitted function with both parameters increased to fit_ap_lp.txt" << std::endl;

    std::ofstream outFilemm("fit_am_lm.txt");
    if (!outFilemm) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return 1;
    };

    for (double t = t_start; t <= t_end; t += step) {
        double value = decay(t, expfit.c[0]-expfit.cerr[0], expfit.c[1]-expfit.cerr[1]);
        outFilemm << t << "\t" << value << "\n";
    };

    outFilemm.close();
    std::cout << "Wrote data for the fitted function with both parameters decreased to fit_am_lm.txt" << std::endl;

    std::ofstream outFilepm("fit_ap_lm.txt");
    if (!outFilepm) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return 1;
    };

    for (double t = t_start; t <= t_end; t += step) {
        double value = decay(t, expfit.c[0]+expfit.cerr[0], expfit.c[1]-expfit.cerr[1]);
        outFilepm << t << "\t" << value << "\n";
    };

    outFilepm.close();
    std::cout << "Wrote data for the fitted function with a increased and lambda decreased to fit_ap_lm.txt" << std::endl;

    std::ofstream outFilemp("fit_am_lp.txt");
    if (!outFilemp) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return 1;
    };

    for (double t = t_start; t <= t_end; t += step) {
        double value = decay(t, expfit.c[0]-expfit.cerr[0], expfit.c[1]+expfit.cerr[1]);
        outFilemp << t << "\t" << value << "\n";
    };

    outFilemp.close();
    std::cout << "Wrote data for the fitted function with a decreased and lambda increased to fit_am_lp.txt" << std::endl;


    return 0;
};



// Test of QR-routine:

// std::random_device rd;
// std::mt19937 gen(rd());
// std::uniform_real_distribution<double> dist(0.1, 10.0);

// int n = 20, m = 8;
// pp::matrix A(n, m); // Tall matrix

// for (int j = 0; j < m; ++j) {
//     std::vector<double> entries;
//     for (int i = 0; i < n; ++i) {
//         entries.push_back(dist(gen));
//     }
//     A.cols[j] = pp::vector(entries);
// }

// std::cout << "Matrix A:" << std::endl;
// A.print();

// QR QR_A(A);
// pp::matrix Q_A = QR_A.Q;
// pp::matrix R_A = QR_A.R;

// std::cout << "Matrix Q:" << std::endl;
// Q_A.print();

// std::cout << "Matrix R:" << std::endl;
// R_A.print();

// std::cout << "Q^T Q:" << std::endl;
// (Q_A.transpose() * Q_A).print();

// std::cout << "Is A=QR? (approximately): " << (A==Q_A*R_A) << std::endl;

// int d = 5;
// pp::matrix B(d, d); // Square matrix
// pp::vector b; // of Ax=b

// for (int j = 0; j < d; ++j) {
//     std::vector<double> entries;
//     for (int i = 0; i < d; ++i) {
//         entries.push_back(dist(gen));
//     }
//     B.cols[j] = pp::vector(entries);
//     b.data.push_back(dist(gen));
// }

// std::cout << "Matrix B:" << std::endl;
// B.print();
// std::cout << "Vector b:" << std::endl;
// b.print();

// QR QR_B(B);
// pp::matrix Q_B = QR_B.Q;
// pp::matrix R_B = QR_B.R;

// std::cout << "Solution x for Bx=b:" << std::endl;
// pp::vector x = QR_B.solve(b);
// x.print();

// std::cout << "Is Bx=b? (approximately): " << (B * x == b) << std::endl;

// std::cout << "Inverse of B:" << std::endl;
// pp::matrix B_inv = QR_B.inverse();
// B_inv.print();

// // Construct the identity matrix
// pp::matrix id_d(d, d);
// for (int i = 0; i < d; i++) {
//     for (int j = 0; j < d; j++) {
//         id_d.set(i, j, i == j ? 1 : 0);
//     }
// }

// std::cout
// << "Is B^(-1) * B and B * B^(-1) the identity matrix? (approximately): "
// << (B_inv * B == id_d) << " and " << (B * B_inv == id_d) << ", respectively."
// << std::endl;