#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <stdexcept> // Used for QR
#include <lapacke.h>
#include <cblas.h> // Used for mat. mult.
#include <fstream> // For writing out the data


bool approx(double a, double b, double abs_prec=1e-9, double rel_prec=1e-9){
    if (std::abs(b-a) <= abs_prec) return true;
    if (std::abs(b-a)/std::max(a, b) <= rel_prec) return true;
    return false;
};

std::vector<double> linspace(double start, double end, int num_points) {
    std::vector<double> result;
    double step = (end - start) / (num_points - 1);  // Calculate step size

    for (int i = 0; i < num_points; ++i) {
        result.push_back(start + i * step);
    }
    return result;
}


// Multiply matrices
std::vector<std::vector<double>> mult_mat(const std::vector<std::vector<double>>& A,
                                                        const std::vector<std::vector<double>>& B) {
    int n = A.size();
    int lda = n;  // Leading dimension for A (number of rows)
    int ldb = n;  // Leading dimension for B (number of rows)
    int ldc = n;  // Leading dimension for C (number of rows)

    // Convert A and B to 1D arrays in column-major order
    std::vector<double> A_flat(n * n);
    std::vector<double> B_flat(n * n);
    std::vector<double> C_flat(n * n, 0.0);  // Initialize result matrix with 0s

    // Copy elements to A_flat in column-major order
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
        A_flat[j * n + i] = A[i][j];
        B_flat[j * n + i] = B[i][j];
    };

    // Perform matrix multiplication C = A * B
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A_flat.data(), lda, B_flat.data(), ldb, 0.0, C_flat.data(), ldc);

    // Convert the result C_flat back to a 2D matrix in row-major order
    std::vector<std::vector<double>> C(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
        C[i][j] = C_flat[j * n + i];
    }

    return C;  // Return the resulting matrix
};

std::vector<std::vector<double>> transposeMatrix(const std::vector<std::vector<double>>& A) {
    int rows = A.size();
    int cols = A[0].size();

    // Create a new matrix of the same size but swapped dimensions (cols x rows)
    std::vector<std::vector<double>> B(cols, std::vector<double>(rows));

    // Transpose the matrix
    for (int i = 0; i < rows; ++i) for (int j = 0; j < cols; ++j) {
        B[j][i] = A[i][j];
    };

    return B;
};


// Jacobi rotation matrix multiplier from the right
void timesJ(std::vector<std::vector<double>>& A, int p, int q, double theta) {
    double c=std::cos(theta), s=std::sin(theta);
    int n = A[0].size(); // Sum over rows of A for matrix mult., just take the first row as an example (should prob. add an error)
    for (int i=0; i < n; i++) {
        double Aip = A[i][p], Aiq = A[i][q];
        A[i][p] = c*Aip - s*Aiq;
        A[i][q] = s*Aip + c*Aiq;
    };
};


// Jacobi rotation matrix multiplier from the left
void Jtimes(std::vector<std::vector<double>>& A, int p, int q, double theta) {
    double c=std::cos(theta), s=std::sin(theta);
    int n = A.size(); // Sum over columns of A for matrix mult.
    for (int j=0; j < n; j++) {
        double Apj = A[p][j], Aqj = A[q][j];
        A[p][j] = c*Apj + s*Aqj;
        A[q][j] = -s*Apj + c*Aqj;
    };
};


// Row-by-row sweeps. Note that this, as well as the mult., change the matrices they're 'fed'
void cycle(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& V) {
    bool updated;
    int n = A.size(); // Only square matrices, so the details of fetching the dim. shouldn't matter
    do {
        updated = false;
        for(int p=0; p<n-1; p++) for(int q=p+1; q < n; q++) { // for(int p=0; p<n; p++) for(int q=0; q < n; q++) if(p != q) {
            double Apq=A[p][q], App=A[p][p], Aqq=A[q][q];
		    double theta=0.5*std::atan2(2*Apq, Aqq-App);
            double c=std::cos(theta), s=std::sin(theta);
            double new_App=c*c*App-2*s*c*Apq+s*s*Aqq;
		    double new_Aqq=s*s*App+2*s*c*Apq+c*c*Aqq;
            if (!approx(App, new_App) || !approx(Aqq, new_Aqq) || !approx(Apq, 0)) { // if (!approx(App, new_App) || !approx(Aqq, new_Aqq)) {
                updated = true;
                timesJ(A, p, q, theta);
                Jtimes(A, p, q, -theta);
                timesJ(V, p, q, theta); // V being the matrix of eigenvectors
            };
        // std::cout << p << q << " " << Apq << std::endl; // Debug print statement
        };
    } while (updated);
    // for(int p=0; p<n-1; p++) for(int q=p+1; q < n; q++) {
    //     A[q][p] = A[p][q]; // Enforce symmetry after "diag."
    // }
};


struct EVD {
    int n; // Only allow square matrices
    std::vector<double> w; // Vector of eigenvalues, since we don't actually construct the diag. matrix
    std::vector<std::vector<double>> A;
    std::vector<std::vector<double>> V;
    std::vector<std::vector<double>> D;
    // Constructor
    EVD(std::vector<std::vector<double>> input_A) { // Should fit with the format D = V^T A V, i.e. A = V D V^T
        n = input_A.size();  // Only allow square matrices
        A = input_A;

        V.resize(n, std::vector<double>(n, 0.0)); // Resize the matrix to nxn, initially filled with 0s
        for (int i = 0; i < n; ++i) {
            V[i][i] = 1.0;
        };

        D = input_A; // Make a copy to diag.
        cycle(D, V);

        w.resize(n, 0.0);
        for (int i = 0; i < n; ++i) {
            w[i] = D[i][i];
        };
    };
};


int main() {
    std::ofstream convfile1("conv_delta.txt");
    if (!convfile1) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    };
    convfile1 << "Δr ε" << std::endl;

    double r_max = 20;
    std::vector<double> delta_rs = linspace(0.05, 1, 1/0.05);
    for (double delta_r : delta_rs) {
        int N = r_max/delta_r - 1; // Matrix size, so careful with the conversion to an integer
        std::vector<double> rs = linspace(delta_r, r_max, N);
        
        // Initialize an NxN matrix
        std::vector<std::vector<double>> H(N, std::vector<double>(N, 0.0));
        
        // Fill the Hamiltonian
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i==j) H[i][j] = 1/(delta_r*delta_r) - 1/rs[i];
                if (i==j+1) H[i][j] = -1/(2*delta_r*delta_r);
                if (i==j-1) H[i][j] = -1/(2*delta_r*delta_r);
            };
        };
        
        EVD solution(H);
        std::vector<double> eigenvals = solution.w;
        double eps = eigenvals[0];

        convfile1 << delta_r << " " << eps << std::endl;
    };
    convfile1.close();
    std::cout << "Delta_r convergence written to conv_delta.txt" << std::endl;


    std::ofstream convfile2("conv_max.txt");
    if (!convfile2) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    };
    convfile2 << "r_max ε" << std::endl;

    double delta_r = 0.05;
    std::vector<double> r_maxs = linspace(5, 30, 30-5+1);
    for (double r_max : r_maxs) {
        int N = r_max/delta_r - 1; // Matrix size, so careful with the conversion to an integer
        std::vector<double> rs = linspace(delta_r, r_max, N);
        
        // Initialize an NxN matrix
        std::vector<std::vector<double>> H(N, std::vector<double>(N, 0.0));
        
        // Fill the Hamiltonian
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i==j) H[i][j] = 1/(delta_r*delta_r) - 1/rs[i];
                if (i==j+1) H[i][j] = -1/(2*delta_r*delta_r);
                if (i==j-1) H[i][j] = -1/(2*delta_r*delta_r);
            };
        };
        
        EVD solution(H);
        std::vector<double> eigenvals = solution.w;
        double eps = eigenvals[0];

        convfile2 << r_max << " " << eps << std::endl;
    };
    convfile2.close();
    std::cout << "r_max convergence written to conv_max.txt" << std::endl;

    return 0;
};