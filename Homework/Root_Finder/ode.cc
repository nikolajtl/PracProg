#include "ode.h"
#include <cmath>
#include <algorithm> // for std::min

namespace ode {

static const double pi = 3.14159265358979323846; // static to limit linkage

std::pair<pp::vector, pp::vector> rkstepXY(
    std::function<pp::vector(double, pp::vector)> f,
    double x,
    pp::vector y,
    double h
) {
    double alpha = 1.0 / 3.0;
    double beta = 2.0 / 3.0;
    pp::vector k0 = f(x, y);
    pp::vector k1 = f(x + alpha * h, y + alpha * k0 * h);
    pp::vector k2 = f(x + beta * h, y + (beta - beta * beta / (2 * alpha)) * k0 * h + beta * beta / (2 * alpha) * k1 * h);

    pp::vector y_next = y + (1 - 1 / (2 * alpha) - 1 / (3 * beta * beta) + 1 / (3 * alpha * beta)) * k0 * h
                        + (1 / (2 * alpha) - 1 / (3 * alpha * beta) + 1 / (3 * beta * beta)) * k1 * h;

    pp::vector y_lower = y + (1 - 1 / (2 * alpha)) * k0 * h + (1 / (2 * alpha)) * k1 * h;
    pp::vector delta_y = (y_next - y_lower) * h;

    return {y_next, delta_y};
}

std::pair<std::vector<double>, std::vector<pp::vector>> driver(
    std::function<pp::vector(double, pp::vector)> f,
    std::pair<double, double> interval,
    pp::vector y_a,
    double h,
    double acc,
    double eps
) {
    std::vector<double> xs;
    std::vector<pp::vector> ys;

    double a = interval.first;
    double b = interval.second;

    xs.push_back(a);
    ys.push_back(y_a);

    double x = a;
    pp::vector y = y_a;

    do {
        if (x >= b) return {xs, ys};
        if (x + h >= b) h = b - x;

        auto [yh, delta_y] = rkstepXY(f, x, y, h);

        double tol = acc + eps * yh.norm() * std::sqrt(h / (b - a));
        double err = delta_y.norm();

        if (err <= tol) {
            x += h;
            y = yh;
            xs.push_back(x);
            ys.push_back(y);
        }
        h *= std::min(std::pow(tol / err, 0.25) * 0.95, 2.0);
    } while (true);
}

} // namespace ode
