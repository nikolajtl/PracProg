#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <stdexcept> // Used for QR
#include <lapacke.h>
#include <cblas.h> // Used for mat. mult.


bool approx(double a, double b, double abs_prec=1e-9, double rel_prec=1e-9){
    if (std::abs(b-a) <= abs_prec) return true;
    if (std::abs(b-a)/std::max(a, b) <= rel_prec) return true;
    return false;
};


// Made with lapacke (by ChatGPT) as to not use my own matrix class from the previous homework
// class QR {
// private:
//     int m, n;                // Dimensions: m rows, n columns
//     std::vector<double> A;   // Original matrix (row-major)
//     std::vector<double> Q;   // Orthogonal matrix
//     std::vector<double> R;   // Upper triangular matrix
//     std::vector<double> tau; // Scalar factors

//     // Accessor for row-major indexing
//     double& access(std::vector<double>& mat, int i, int j, int cols) {
//         return mat[i * cols + j];
//     }

// public:
//     QR(const std::vector<std::vector<double>>& inputA) : m(inputA.size()), n(inputA[0].size()), A(m * n), Q(m * n), R(n * n, 0), tau(n) {
//         // Flatten input matrix (row-major order)
//         for (int i = 0; i < m; ++i) {
//             for (int j = 0; j < n; ++j) {
//                 A[i * n + j] = inputA[i][j];
//             }
//         }

//         // Transpose A to column-major for LAPACKE
//         std::vector<double> A_col_major(m * n);
//         for (int i = 0; i < m; ++i) {
//             for (int j = 0; j < n; ++j) {
//                 A_col_major[j * m + i] = A[i * n + j];
//             }
//         }

//         // Compute QR decomposition
//         int info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, m, n, A_col_major.data(), m, tau.data());
//         if (info != 0) throw std::runtime_error("QR decomposition failed");

//         // Extract R (upper triangular) from A_col_major
//         for (int i = 0; i < n; ++i) {
//             for (int j = i; j < n; ++j) {
//                 access(R, i, j, n) = A_col_major[j * m + i];
//             }
//         }

//         // Generate orthogonal matrix Q
//         info = LAPACKE_dorgqr(LAPACK_COL_MAJOR, m, n, n, A_col_major.data(), m, tau.data());
//         if (info != 0) throw std::runtime_error("Failed to compute Q matrix");

//         // Store Q in row-major format
//         for (int i = 0; i < m; ++i) {
//             for (int j = 0; j < n; ++j) {
//                 access(Q, i, j, n) = A_col_major[j * m + i];
//             }
//         }
//     }

//     std::vector<double> solve(const std::vector<double>& b) {
//         if (b.size() != static_cast<size_t>(m)) throw std::runtime_error("Dimension mismatch");

//         std::vector<double> x(b);

//         // Transform b to column-major for LAPACKE
//         std::vector<double> b_col_major = x;

//         // Compute Q^T * b
//         int info = LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'T', m, 1, n, Q.data(), m, tau.data(), b_col_major.data(), m);
//         if (info != 0) throw std::runtime_error("Failed to compute Q^T * b");

//         // Solve Rx = Q^T * b
//         info = LAPACKE_dtrtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', n, 1, R.data(), n, b_col_major.data(), n);
//         if (info != 0) throw std::runtime_error("Solving Rx = c failed");

//         return std::vector<double>(b_col_major.begin(), b_col_major.begin() + n);
//     }

//     double det() {
//         double determinant = 1.0;
//         for (int i = 0; i < n; ++i) {
//             determinant *= access(R, i, i, n);
//         }
//         return determinant;
//     }

//     std::vector<std::vector<double>> inverse() {
//         if (m != n) throw std::runtime_error("Matrix must be square");

//         std::vector<std::vector<double>> inv(n, std::vector<double>(n, 0));

//         // Solve A * X = I column by column
//         for (int j = 0; j < n; ++j) {
//             std::vector<double> e(n, 0);
//             e[j] = 1.0;
//             std::vector<double> x = solve(e);

//             for (int i = 0; i < n; ++i) {
//                 inv[i][j] = x[i];
//             }
//         }

//         return inv;
//     }

//     void print_matrix(const std::vector<double>& mat, int rows, int cols) {
//         for (int i = 0; i < rows; ++i) {
//             for (int j = 0; j < cols; ++j) {
//                 std::cout << mat[i * cols + j] << " ";
//             }
//             std::cout << "\n";
//         }
//     }
// };

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
        A[q][j] = s*Apj + c*Aqj;
    };
};


// Row-by-row sweeps. Note that this, as well as the mult., change the matrices they're 'fed'
void cycle(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& V) {
    bool updated;
    do {
        updated = false;
        int n = A.size(); // Only square matrices, so the details of fetching the dim. shouldn't matter
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
        std::cout << p << q << " " << Apq << std::endl;
        };
    } while (updated);
};


struct EVD {
    int n; // Only allow square matrices
    std::vector<double> w; // Vector of eigenvalues, IDK why one should make one-such
    std::vector<std::vector<double>> A;
    std::vector<std::vector<double>> V;
    std::vector<std::vector<double>> D;
    // Constructor
    EVD(std::vector<std::vector<double>> input_A) { // Should fit with the format D = V^T A V, i.e. A = V D V^T
        n = A.size();  // Only allow square matrices
        A = input_A;

        V.resize(n, std::vector<double>(n, 0.0)); // Resize the matrix to nxn, initially filled with 0s
        for (int i = 0; i < n; ++i) {
            V[i][i] = 1.0;
        };

        std::vector<std::vector<double>> D = input_A; // Make a copy to diag.
        cycle(D, V);

        w.resize(n, 0.0);
        for (int i = 0; i < n; ++i) {
            w[i] = D[i][i];
        };
    };
};


int main(int argc, char* argv[]) {
    // Ensure proper argument usage
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <size of matrix>" << std::endl;
        return 1;
    }

    int N = std::atoi(argv[1]); // Get matrix size from command line argument
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.1, 10.0); // Random values between 0.1 and 10.0

    // Initialize an NxN matrix
    std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0));

    // Fill the matrix to be symmetric
    for (int i = 0; i < N; ++i) {
        for (int j = i; j < N; ++j) {
            double value = dist(gen);
            A[i][j] = value;     // Set the upper triangular part
            A[j][i] = value;     // Set the lower triangular part (symmetry)
        };
    };

    // Print the matrix to verify it's symmetric
    std::cout << "Generated symmetric matrix A:" << std::endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    };

    // EVD EVD_A(A);
    // std::vector<std::vector<double>> D = EVD_A.D;
    // std::vector<std::vector<double>> V = EVD_A.V;
    // std::vector<std::vector<double>> V_T = transposeMatrix(V);
    // std::vector<std::vector<double>> VTAV = mult_mat(mult_mat(V_T, A), V);
    // std::vector<std::vector<double>> VDVT = mult_mat(mult_mat(V, D), V_T);
    // std::vector<std::vector<double>> VTV = mult_mat(V_T, V);
    // std::vector<std::vector<double>> VVT = mult_mat(V, V_T);

    // std::cout << "V^T V and V V^T, respectively. Should be two identities:" << std::endl;
    // for (int i = 0; i < N; ++i) {
    //     for (int j = 0; j < N; ++j) {
    //         std::cout << VTV[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // };
    // for (int i = 0; i < N; ++i) {
    //     for (int j = 0; j < N; ++j) {
    //         std::cout << VVT[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // };

    std::vector<std::vector<double>> V;
    V.resize(N, std::vector<double>(N, 0.0)); // Resize the matrix to nxn, initially filled with 0s
        for (int i = 0; i < N; ++i) {
            V[i][i] = 1.0;
        };
    
    std::cout << "V before cycle:" << std::endl;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                std::cout << V[i][j] << " ";
            }
            std::cout << std::endl;
    };

    cycle(A, V);

    std::cout << "A after cycle:" << std::endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    };
    
    std::cout << "V after cycle:" << std::endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << V[i][j] << " ";
        }
        std::cout << std::endl;
    };

    std::vector<std::vector<double>> V_T = transposeMatrix(V);
    std::vector<std::vector<double>> VTV = mult_mat(V_T, V);
    std::cout << "V^T V:" << std::endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << VTV[i][j] << " ";
        }
        std::cout << std::endl;
    };

    return 0;
};


// Eksempel på brug:
// int main() {
//     std::vector<std::vector<double>> A = {
//         {12, -51, 4},
//         {6, 167, -68},
//         {-4, 24, -41}
//     };

//     return 0;
// };