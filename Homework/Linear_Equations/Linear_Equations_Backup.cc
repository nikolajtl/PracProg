#include <iostream>
#include <vector>
#include <cmath>
#include <random> // Random numbers

struct vect{
    std::vector<double> entries; // Entries of the vector
    int dim; // Dimension of the vector space

    // Constructor
    vect() {};
    vect(std::vector<double> entries) : entries(entries), dim(entries.size()) {};

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

    double dot(const vect& other){
        double dot_product = 0; // Initialize the dot product
        for(int i = 0; i < this->dim; i++){
            dot_product += this->entries[i]*other.entries[i]; // Add the product of the entries
        }
        return dot_product; // Return the dot product
    }

    double norm(){
        double normsquare = dot(*this);
        double norm = sqrt(normsquare);
        return norm;
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
};


struct matrix{
    std::vector<vect> columns; // Define the matrix as a set of column vectors
    vect diag; // Define the diagonal of the matrix
    int colnr; // Dimension of the column space
    int rownr; // Dimension of the row space

    // Constructors
    matrix() {};
    matrix(std::vector<vect> columns) : columns(columns), colnr(columns.size()), rownr(columns[0].dim) {
        for (const auto& col : columns){
            if (col.dim != rownr) {
                throw std::invalid_argument("All columns must have the same dimension");
            }
        }
    }
    matrix(vect diag) : diag(diag), colnr(diag.dim), rownr(diag.dim) {
        for (int j=0; j < diag.dim; j++) {
            vect column;
            double d = diag.get(j);
            for (int i = 0; i < diag.dim; i++) {
                if (i == j) {
                    column.entries.push_back(d);
                } else {
                    column.entries.push_back(0);
                }
            }
            columns.push_back(column);
        }
    }

    // Getter: Get the element at row 'i' and column 'j'
    double get(int i, int j) const {
        if (i < 0 || i >= rownr || j < 0 || j >= colnr) {
            throw std::out_of_range("Index out of range");
        }
        return columns[j].entries[i]; // Return the element at the correct position in the column
    }

    // Setter: Set the element at row 'i' and column 'j'
    void set(int i, int j, double value) {
        if (i < 0 || i >= rownr || j < 0 || j >= colnr) {
            throw std::out_of_range("Index out of range");
        }
        columns[j].entries[i] = value; // Set the element at the correct position in the column
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
        std::vector<vect> new_columns;
        for (int i = 0; i < other.colnr; i++) {
            vect new_column = *this * other.columns[i];
            new_columns.push_back(new_column);
        }
        return matrix(new_columns);
    }

    // Transpose
    matrix transpose() {
        std::vector<vect> new_columns;
        for (int i = 0; i < rownr; i++) {
            std::vector<double> new_entries;
            for (int j = 0; j < colnr; j++) {
                new_entries.push_back(columns[j].entries[i]);
            }
            new_columns.push_back(vect(new_entries));
        }
        columns = new_columns;
        return matrix(new_columns);
    }

    // Print to console
    void print() {
        for (int i = 0; i < rownr; i++) {
            for (int j = 0; j < colnr; j++) {
                std::cout << get(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }
};


struct QR{
    matrix A;
    matrix Q;
    matrix R;

    int n; // Number of columns of A (to be)

    QR(matrix A) : A(A), n(A.colnr) { // Note the lack of fail-safe for lin. dependent columns
        R.rownr = A.colnr;
        R.colnr = A.colnr;

        Q.columns.clear();
        R.columns.clear();

        Q.rownr = A.rownr;
        Q.colnr = n;

        for (int i=0; i < n; i++) {
            vect a = A.columns[i];
            vect q = a / a.norm();
            Q.columns.push_back(q);

            std::vector<double> r_entries;
            for (int j = 0; j < n; j++) {
                if (j >= i) { // Combined upper triangular and j = i + 1 to n
                    double projection = A.columns[j].dot(q); // The inner product of a_j with q_i
                    r_entries.push_back(projection);
                    A.columns[j] = A.columns[j] - q * projection;
                } else {
                    r_entries.push_back(0); // Fill with zeros for elements above the diagonal
                }
            }
            R.columns.push_back(vect(r_entries));
        }
    }

    vect solve(vect b) { // Solve the system Ax = b by back substitution on Rx=Q^T*b
        vect x, x_opp;
        matrix U = R;
        vect c = Q.transpose() * b;
        for (int i = c.dim; i >= 0; i--) {
            double U_ii = U.get(i, i);
            double c_i = c.get(i);
            double sum = 0;
            for (int k=i+1; k <= c.dim; k++) {
                sum += U.get(i, k) * x.get(k);
            }
            double x_i = (c_i - sum) / U_ii;
            x_opp.entries.push_back(x_i);
        }
        for (int i = x_opp.dim; i > 0; i--) {
            x.entries.push_back(x_opp.get(i));
        }
        return x;
    }

    double det() {
        double det = 1;
        for (int i = 0; i < n; i++) {
            det *= R.get(i, i);
        }
        return det;
    }
};


int main(){
    std::random_device rd; 
    std::mt19937 gen(rd()); // Use Mersenne Twister engine
    std::uniform_real_distribution<double> dist(0.1, 10.0); // Define a distribution between 0.1 and 10

    int n = 20; // Number of columns; the dimension of the column vectors
    int m = 8; // Number of rows, chosen s.t. n > m
    vect a;
    matrix A;
    for (int i=0; i < m; i++){
        for (int j=0; j < n; j++){
            double val = dist(gen);
            a.entries.push_back(val);
        }
        A.columns.push_back(a);
    }

    QR Q_R(A);
    matrix Q = Q_R.Q;
    matrix R = Q_R.R;
    R.print();
    std::cout << Q.rownr << std::endl;

    return 0;
};