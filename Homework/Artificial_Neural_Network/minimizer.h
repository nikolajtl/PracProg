#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <functional>
#include <stdexcept>
#include "matrix.h"

namespace min {

pp::vector gradient(std::function<double(pp::vector)> φ, pp::vector x, int mode = 0);
pp::matrix hessian(std::function<double(pp::vector)> φ, pp::vector x, int mode = 0);
std::pair<pp::vector, int> newton(std::function<double(pp::vector)> φ, pp::vector x, double acc = 1e-3, int mode = 0);

struct QR {
    pp::matrix A, Q, R;

    QR(const pp::matrix& A);

    pp::vector solve(const pp::vector& b);
    double det();
    pp::matrix inverse();
    pp::matrix Rinverse();
};

} // namespace min

#endif // MINIMIZER_H
