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



// Jacobi rotation matrix multiplier from the right
void timesJ(pp::matrix& A, int p, int q, double theta) {
    double c=std::cos(theta), s=std::sin(theta);
    int n = A[0].size(); // Sum over rows of A for matrix mult., just take the first row as an example (should prob. add an error)
    for (int i=0; i < n; i++) {
        double Aip = A[i][p], Aiq = A[i][q];
        A[i][p] = c*Aip - s*Aiq;
        A[i][q] = s*Aip + c*Aiq;
    };
};


// Jacobi rotation matrix multiplier from the left
void Jtimes(pp::matrix& A, int p, int q, double theta) {
    double c=std::cos(theta), s=std::sin(theta);
    int n = A.size1(); // Sum over columns of A for matrix mult.
    for (int j=0; j < n; j++) {
        double Apj = A[p][j], Aqj = A[q][j];
        A[p][j] = c*Apj + s*Aqj;
        A[q][j] = -s*Apj + c*Aqj;
    };
};


struct EVD {
    int n; // Only allow square matrices
    pp::vector w; // Vector of eigenvalues, since we don't actually construct the diag. matrix
    pp::matrix A;
    pp::matrix V;
    pp::matrix D;
    // Constructor
    EVD(pp::matrix input_A) { // Should fit with the format D = V^T A V, i.e. A = V D V^T
        n = input_A.size1();  // Only allow square matrices
        A = input_A; // Make a copy since the mult. functions pass by reference, preserving A here and making a copy D later to modify

        V.resize(n, n); // Resize the matrix to nxn, initially filled with 0s
        for (int i = 0; i < n; ++i) {
            V[i][i] = 1.0;
        };

        D = input_A; // Make a copy to diagonalize
        int max_sweeps = 5000;
        double acc = 1e-12;
        
        for (int sweep = 0; sweep < max_sweeps; sweep++) {
            bool converged = true;

            for (int p=0; p<n-1; p++) {
                int q = p + 1;
                
                double App = D[p][p], Aqq = D[q][q], Apq = D[p][q];
                
                if (std::abs(Apq) < acc)
                continue;
                
                converged = false;
                
                double theta=0.5*std::atan2(2*Apq, Aqq-App);
                
                timesJ(D, p, q, theta);
                Jtimes(D, p, q, -theta);
                
                timesJ(V, p, q, theta); // Accumulate eigenvectors
            };
            
            if (converged) break;
        };


        w.resize(n);
        for (int i = 0; i < n; ++i) {
            w[i] = D[i][i];
        };
    };
};


int main() {
    int N = 10;

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

    H.print();

    EVD decomp(H);
    pp::vector eigenvals = decomp.w;
    // std::cout << "Found eigenvalues (energies) for the Hamiltonian: " << std::endl;
    // eigenvals.print();
    for (int i=0; i<10; i++) std::cout << eigenvals[i] << std::endl;
    decomp.D.print();

    return 0;
};