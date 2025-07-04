#include <iostream>
#include <fstream>
#include <cmath>
#include <utility>
#include <random>
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
    
    // Debugging
    // std::cout << std::endl;
    // pp::vector test(alpha);
    // test.print();

    return std::make_pair(T, Q);
};



int main() {
    int N = 8; // A is an N x N symmetric matrix
    int n = 6; // n =< N, the resulting T_n is an n x n matrix

    // Make a random q_1 and symmetric A:
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.1, 10.0); // Random values between 0.1 and 10.0
    
    pp::vector q_1(N);
    for (int i=0; i<q_1.size(); i++) q_1[i] = dist(gen);
    std::cout << "Some random q_1: " << std::endl;
    q_1.print();

    pp::matrix A(N, N);
    for (int i=0; i<N; i++) {
        for (int j=i; j<N; j++) {
            double value = dist(gen);
            A[i][j] = value;     // Set the upper triangular part
            A[j][i] = value;     // Set the lower triangular part (symmetry)
        };
    };

    std::cout << "Some random symmetric N x N matrix A for N=" << N << std::endl;
    A.print();

    auto decomp = lanczos(A, q_1, n);
    std::cout << "The Lanczos decomposition T of A for n=" << n << std::endl;
    decomp.first.print();

    pp::matrix Qtest = decomp.second;

    std::cout << "Computed Q^T A Q: " << std::endl;
    pp::matrix Ttest = Qtest.transpose() * A * Qtest;
    Ttest.print();

    std::cout << "Computed Q^T Q for sanity check (should be identity): " << std::endl;
    pp::matrix Itest = Qtest.transpose() * Qtest;
    Itest.print();



    return 0;
};