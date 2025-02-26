#include <iostream>
#include <vector>
#include <cmath>
#include <random>

bool approx(double a, double b, double abs_prec=1e-9, double rel_prec=1e-9){
    if (std::abs(b-a) <= abs_prec) return true;
    if (std::abs(b-a)/std::max(a, b) <= rel_prec) return true;
    return false;
}

struct vect {
    std::vector<double> entries;
    int dim;

    vect() : dim(0) {}
    vect(std::vector<double> entries) : entries(entries), dim(entries.size()) {}

    // Dårlig løsning:
    void update_dim() {
        dim = entries.size();
    }

    double get(int i) const {
        if (i < 0 || i >= dim) {
            throw std::out_of_range("Index out of range");
        }
        return entries[i];
    }

    void set(int i, double value) {
        if (i < 0 || i >= dim) {
            throw std::out_of_range("Index out of range");
        }
        entries[i] = value;
    }

    double dot(const vect& other) const {
        double result = 0;
        for (int i = 0; i < dim; ++i) {
            result += entries[i] * other.entries[i];
        }
        return result;
    }

    double norm() const {
        return std::sqrt(dot(*this));
    }

    vect operator+(const vect& other) const {
        if (this->dim != other.dim) {
            throw std::invalid_argument("Vectors must be of the same dimension");
        }
        std::vector<double> result_entries;
        for (int i = 0; i < this->dim; ++i) {
            result_entries.push_back(this->entries[i] + other.entries[i]);
        }
        return vect(result_entries);
    }

    vect operator-(const vect& other) const {
        if (this->dim != other.dim) {
            throw std::invalid_argument("Vectors must be of the same dimension");
        }
        std::vector<double> result_entries;
        for (int i = 0; i < this->dim; ++i) {
            result_entries.push_back(this->entries[i] - other.entries[i]);
        }
        return vect(result_entries);
    }

    vect operator*(double scalar) const {
        std::vector<double> result_entries;
        for (int i = 0; i < this->dim; ++i) {
            result_entries.push_back(this->entries[i] * scalar);
        }
        return vect(result_entries);
    }

    vect operator/(double scalar) const {
        if (scalar == 0) {
            throw std::invalid_argument("Cannot divide by zero");
        }
        std::vector<double> result_entries;
        for (int i = 0; i < this->dim; ++i) {
            result_entries.push_back(this->entries[i] / scalar);
        }
        return vect(result_entries);
    }

    // Left scalar multiplication
    friend vect operator*(double scalar, const vect& vec) {
        return vec * scalar;
    }

    // Compare two vectors entry-wise
    bool operator==(const vect& other) const {
        if (dim != other.dim) {
            return false;
        }
        for (int i = 0; i < dim; ++i) {
            if (!approx(entries[i], other.entries[i])) {
                return false;
            }
        }
        return true;
    }

    void print() const {
        for (double val : entries) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
};

struct matrix {
    std::vector<vect> columns;
    int rownr, colnr;
    vect diag;

    matrix(int rows, int cols) : rownr(rows), colnr(cols) {
        columns.resize(cols, vect(std::vector<double>(rows, 0)));
    }

    double get(int i, int j) const {
        return columns[j].entries[i];
    }

    void set(int i, int j, double value) {
        columns[j].entries[i] = value;
    }

    // Matrix-vector multiplication
    vect operator*(const vect& v) const {
        if (v.dim != colnr) {
            throw std::invalid_argument("Matrix and vector dimensions do not match");
        }
        std::vector<double> result_entries;
        for (int i = 0; i < rownr; i++) {
            double result = 0;
            for (int j = 0; j < colnr; ++j) {
                result += columns[j].entries[i] * v.entries[j];
            }
            result_entries.push_back(result);
        }
        return vect(result_entries);
    }

    // Matrix-matrix multiplication
    matrix operator*(const matrix& other) const {
        if (colnr != other.rownr) {
            throw std::invalid_argument("Matrix dimensions do not match");
        }
        matrix new_matrix(rownr, other.colnr);
        for (int i = 0; i < other.colnr; i++) {
            vect new_column = *this * other.columns[i];
            new_matrix.columns[i] = new_column;
        }
        return new_matrix;
    }

    // Matrix-matrix subtraction
    matrix operator-(const matrix& other) const {
        if (rownr != other.rownr || colnr != other.colnr) {
            throw std::invalid_argument("Matrix dimensions do not match");
        }
        matrix new_matrix(rownr, colnr);
        for (int i = 0; i < rownr; i++) {
            for (int j = 0; j < colnr; j++) {
                new_matrix.set(i, j, get(i, j) - other.get(i, j));
            }
        }
        return new_matrix;
    }

    // Transpose
    matrix transpose() const {
        matrix new_matrix(colnr, rownr);
        for (int i = 0; i < rownr; i++) {
            for (int j = 0; j < colnr; j++) {
                new_matrix.set(j, i, get(i, j));
            }
        }
        return new_matrix;
    }

    // Print to console
    void print() const {
        for (int i = 0; i < rownr; ++i) {
            for (int j = 0; j < colnr; ++j) {
                std::cout << get(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }

    // Compare two matrices entry-wise
    bool operator==(const matrix& other) const {
        if (rownr != other.rownr || colnr != other.colnr) {
            return false;
        }
        for (int i = 0; i < rownr; ++i) {
            for (int j = 0; j < colnr; ++j) {
                if (!approx(get(i, j), other.get(i, j))) {
                    return false;
                }
            }
        }
        return true;
    }
};

struct QR {
    matrix A, Q, R;

    QR(const matrix& A) : A(A), Q(A.rownr, A.colnr), R(A.colnr, A.colnr) {
        for (int i = 0; i < A.colnr; ++i) {
            vect a = A.columns[i];

            for (int j = 0; j < i; ++j) {
                double r = a.dot(Q.columns[j]);
                R.set(j, i, r);
                a = a - (Q.columns[j] * r);
            }

            double norm = a.norm();
            if (norm == 0) throw std::runtime_error("Linear dependency detected");

            Q.columns[i] = a / norm;
            R.set(i, i, norm);
        }
    }

    vect solve(const vect& b) {
        vect c = Q.transpose() * b;

        std::vector<double> x_vals(R.rownr, 0);
        for (int i = R.rownr - 1; i >= 0; --i) {
            double sum = 0;
            for (int j = i + 1; j < R.rownr; ++j) {
                sum += R.get(i, j) * x_vals[j];
            }
            x_vals[i] = (c.entries[i] - sum) / R.get(i, i);
        }
        return vect(x_vals);
    }

    double det() {
        int n = R.rownr;
        double det = 1;
        for (int i = 0; i < n; i++) {
            det *= R.get(i, i);
        }
        return det;
    }

    matrix inverse() {
        if (det() == 0) throw std::runtime_error("Matrix is singular");
        if (A.rownr != A.colnr) throw std::runtime_error("Matrix is not square");

        matrix inv(A.colnr, A.rownr);
        for (int i = 0; i < A.colnr; i++) {
            // Generate the i-th unit vector
            vect e_i;
            for (int j = 0; j < A.rownr; j++) {
                e_i.entries.push_back(i == j ? 1 : 0);
            }
            e_i.update_dim();

            // Solve A * x_i = e_i
            vect x = solve(e_i);
            for (int j = 0; j < A.rownr; j++) {
                inv.set(j, i, x.get(j)); // x_i is the i'th column of the inverse
            }
        }
        return inv;
    }
};

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.1, 10.0);

    int n = 20, m = 8;
    matrix A(n, m); // Tall matrix

    for (int j = 0; j < m; ++j) {
        std::vector<double> entries;
        for (int i = 0; i < n; ++i) {
            entries.push_back(dist(gen));
        }
        A.columns[j] = vect(entries);
    }

    std::cout << "Matrix A:" << std::endl;
    A.print();

    QR QR_A(A);
    matrix Q_A = QR_A.Q;
    matrix R_A = QR_A.R;

    std::cout << "Matrix Q:" << std::endl;
    Q_A.print();

    std::cout << "Matrix R:" << std::endl;
    R_A.print();

    std::cout << "Q^T Q:" << std::endl;
    (Q_A.transpose() * Q_A).print();

    std::cout << "Is A=QR? (approximately): " << (A==Q_A*R_A) << std::endl;

    int d = 5;
    matrix B(d, d); // Square matrix
    vect b; // of Ax=b

    for (int j = 0; j < d; ++j) {
        std::vector<double> entries;
        for (int i = 0; i < d; ++i) {
            entries.push_back(dist(gen));
        }
        B.columns[j] = vect(entries);
        b.entries.push_back(dist(gen));
    }
    b.update_dim();

    std::cout << "Matrix B:" << std::endl;
    B.print();
    std::cout << "Vector b:" << std::endl;
    b.print();

    QR QR_B(B);
    matrix Q_B = QR_B.Q;
    matrix R_B = QR_B.R;

    std::cout << "Solution x for Bx=b:" << std::endl;
    vect x = QR_B.solve(b);
    x.print();

    std::cout << "Is Bx=b? (approximately): " << (B * x == b) << std::endl;

    std::cout << "Inverse of B:" << std::endl;
    matrix B_inv = QR_B.inverse();
    B_inv.print();

    // Construct the identity matrix
    matrix id_d(d, d);
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            id_d.set(i, j, i == j ? 1 : 0);
        }
    }

    std::cout
    << "Is B^(-1) * B and B * B^(-1) the identity matrix? (approximately): "
    << (B_inv * B == id_d) << " and " << (B * B_inv == id_d) << ", respectively."
    << std::endl;

    return 0;
}