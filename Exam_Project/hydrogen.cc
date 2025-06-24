#include <iostream>
#include <fstream>
#include <cmath>
#include <utility>
#include <random>
#include "matrix.h"



std::vector<double> linspace(double start, double end, int num_points) {
    std::vector<double> result;
    double step = (end - start) / (num_points - 1);  // Calculate step size

    for (int i = 0; i < num_points; ++i) {
        result.push_back(start + i * step);
    }
    return result;
};



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



// Takes matrix A=A_k as input, decomposes A_k = Q_k R_k, and returns A_{k+1} = R_k Q_k
pp::matrix QR_updater(pp::matrix A) {
    QR decomp(A);
    pp::matrix Q = decomp.Q;
    pp::matrix R = decomp.R;

    return (R * Q);
};

// Uses the updater for at most k_max iterations to diagonalize the input A, given some precision
pp::matrix QR_diag(pp::matrix A, double acc, int k_max=1000) {
    int Adim = A.size1();
    if (Adim != A.size2()) std::cerr << "Error: Attempted diagonalization of non-square matrix" << std::endl;

    pp::matrix A_k = A;
    for (int k=0; k<k_max; k++) {
        A_k = QR_updater(A_k);

        // Check convergence: are the off-diagonal elements small?
        bool converged = true;
        for (int i=0; i < Adim-1; i++) {
            if (std::abs(A_k(i + 1, i)) > acc) {
                converged = false;
                break;
            }
        }
        
        if (converged) break;

        // Debugging
        // std::cout << k << std::endl;
    };

    return A_k;
};



// First entry T, second Q
std::pair<pp::matrix, pp::matrix> lanczos(const pp::matrix& A, pp::vector q_1, int dim) {
    int n = A.size1();
    q_1 /= q_1.norm(); // Normalize the starting vector

    pp::matrix Q(n, dim);
    Q[0] = q_1;
    pp::vector q_k = q_1;
    pp::vector q_prev = q_1; // For sending along to the loop

    std::vector<NUMBER> alpha;
    std::vector<NUMBER> beta;
    double alpha_k, beta_k;

    // Iterative steps (k = 2 to dim but zero-indexed)
    for (int k = 1; k < dim; k++) {
        q_k = A * q_prev;
        alpha_k = q_k.dot(q_prev);

        for (int i=0; i<k; i++) {
            pp::vector v = Q[i];
            double coeff = q_k.dot(v);
            q_k -= coeff * v;
        };

        beta_k = q_k.norm();

        alpha.push_back(alpha_k);
        beta.push_back(beta_k);

        q_k /= beta_k;
        Q[k] = q_k;

        q_prev = q_k;
    };

    // Last step for the lower-right (n, n) diagonal (alpha) entry:
    q_k = A* q_prev;
    alpha_k = q_k.dot(q_prev);
    alpha.push_back(alpha_k);

    // Build the tridiagonal matrix T
    pp::matrix T(dim, dim);
    for (int i=0; i<dim; i++) {
        T(i, i) = alpha[i];
        if (i + 1 < dim) {
            T(i, i + 1) = beta[i];
            T(i + 1, i) = beta[i];
        };
    };

    return std::make_pair(T, Q);
};



int main() {
    int N = 200; // A (here H) is an N x N symmetric matrix
    int n = 10; // n =< N, the resulting T_n is an n x n matrix

    double delta_r = 0.1;
    double r_max = (N+1)*delta_r; // N+1 to fit with the boundary conditions as outlined in the eigenvalue homework
    std::vector<double> rs = linspace(delta_r, r_max, N);

    // Make the Hamiltonian
    pp::matrix H(N, N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i==j) H[i][j] = 1/(delta_r*delta_r) - 1/rs[i];
            if (i==j+1) H[i][j] = -1/(2*delta_r*delta_r);
            if (i==j-1) H[i][j] = -1/(2*delta_r*delta_r);
        };
    };
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.1, 1.0); // Random values between 0.1 and 1.0

    // q_1 is here taken to be the analytic ground state with random added noise
    pp::vector q_1(N);
    for (int i=0; i<q_1.size(); i++) q_1[i] = dist(gen) * rs[i] * std::exp(-rs[i]);

    auto decomp = lanczos(H, q_1, n);
    pp::matrix T = decomp.first;
    std::cout << "The n x n for n=" << n << " Lanczos T matrix of our N x N for N=" << N << " Hamiltonian: " << std::endl;
    T.print();

    pp::matrix T_diag = QR_diag(T, 1e-12);
    std::cout << "The QR-diagonalized T: " << std::endl;
    T_diag.print();

    // Convergence check
    std::vector<int> ns;
    for (int i=5; i<20; i++) ns.push_back(i);

    std::ofstream convfile("convergence.txt");
    for (int n_loop : ns) {
        auto decomp_loop = lanczos(H, q_1, n_loop);
        pp::matrix T_n = decomp_loop.first;
        pp::matrix T_n_diag = QR_diag(T_n, 1e-12);
        convfile << n_loop << " " << T_n_diag[n_loop-1][n_loop-1] << std::endl;
    };
    convfile.close();

    return 0;
};