#include <cmath>
#include <stdexcept>
#include "minimizer.h"
#include "matrix.h"

#include <iostream>

namespace min {

QR::QR(const pp::matrix& A) : A(A), Q(A.size1(), A.size2()), R(A.size2(), A.size2()) {
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

pp::vector QR::solve(const pp::vector& b) {
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

double QR::det() {
    int n = R.size1();
    double det = 1;
    for (int i = 0; i < n; i++) {
        det *= R.get(i, i);
    }
    return det;
};

pp::matrix QR::inverse() {
    if (det() == 0) throw std::runtime_error("Matrix is singular");
    if (A.size1() != A.size2()) throw std::runtime_error("Matrix is not square");

    pp::matrix inv(A.size2(), A.size1());
    for (int i = 0; i < A.size2(); i++) {
        pp::vector e_i;
        for (int j = 0; j < A.size1(); j++) {
            e_i.data.push_back(i == j ? 1 : 0);
        }

        pp::vector x = solve(e_i);
        for (int j = 0; j < A.size1(); j++) {
            inv.set(j, i, x[j]);
        }
    }
    return inv;
};

pp::matrix QR::Rinverse() {
    if (R.size1() != R.size2()) throw std::runtime_error("R is not square");

    int n = R.size1();
    pp::matrix Rinv(n, n);

    for (int i = n - 1; i >= 0; --i) {
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

pp::vector gradient(std::function<double(pp::vector)> φ, pp::vector x, int mode) {
    double φx = φ(x), φx_plus = φx, φx_minus = φx;
    pp::vector gφ_forward(x.size()), gφ_central(x.size());

    for (int i = 0; i < x.size(); ++i) {
        // double dxi = std::abs(x[i]) * std::pow(2, -26);
        double dxi = std::abs(x[i]) > 1e-8 ? std::abs(x[i]) * std::pow(2, -26) : 1e-8;
        x[i] += dxi;
        φx_plus = φ(x);
        gφ_forward[i] = (φx_plus - φx) / dxi;

        x[i] -= 2 * dxi;
        φx_minus = φ(x);
        gφ_central[i] = (φx_plus - φx_minus) / (2 * dxi);

        x[i] += dxi; // Restore original
    }

    if (mode == 0)
        return gφ_forward;
    else if (mode == 1)
        return gφ_central;
    else
        throw std::runtime_error("Invalid mode in gradient; must be 0 (forward) or 1 (central)");
}


pp::matrix hessian(std::function<double(pp::vector)> φ, pp::vector x, int mode) {
    pp::matrix H_forward(x.size(), x.size());
    pp::matrix H_central(x.size(), x.size());
    pp::vector gφx_forward = gradient(φ, x, 0);
    pp::vector gφx_central = gradient(φ, x, 1);
    pp::vector gφx_forward_plus, gφx_central_plus, gφx_central_minus;

    for (int j = 0; j < x.size(); ++j) {
        // double dxj = std::abs(x[j]) * std::pow(2, -26);
        double dxj = std::abs(x[j]) > 1e-8 ? std::abs(x[j]) * std::pow(2, -26) : 1e-8;

        x[j] += dxj;
        gφx_forward_plus = gradient(φ, x, 0);
        gφx_central_plus = gradient(φ, x, 1);

        x[j] -= 2 * dxj;
        gφx_central_minus = gradient(φ, x, 1);

        x[j] += dxj;

        for (int i = 0; i < x.size(); ++i) {
            H_forward(i, j) = (gφx_forward_plus[i] - gφx_forward[i]) / dxj;
            H_central(i, j) = (gφx_central_plus[i] - gφx_central_minus[i]) / (2 * dxj);
        }
    }

    if (mode == 0)
        return H_forward;
    else if (mode == 1)
        return H_central;
    else
        throw std::runtime_error("Invalid mode in hessian; must be 0 (forward) or 1 (central)");
}


std::pair<pp::vector, int> newton(std::function<double(pp::vector)> φ, pp::vector x, double acc, int mode) {
    int steps = 0;

    while (true) {
        pp::vector g = gradient(φ, x, mode);

        // Debugging
        // std::cout << g.norm() << std::endl;
        std::cout << φ(x) << std::endl;

        if (φ(x) < acc)
            break;

        pp::matrix H = hessian(φ, x, mode);
        QR qrH(H);
        pp::vector dx = qrH.solve(-1 * g);

        double λ = 1.0;
        while (λ > 1.0 / 128.0) {
            if (φ(x + λ * dx) < φ(x))
                break;
            λ /= 2;
        }

        x = x + λ * dx;
        ++steps;

        // if (steps > 1e3) break;
    }

    return {x, steps};
}

} // namespace min