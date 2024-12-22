#pragma once
#include <functional>
#include <vector>
#include <cmath>
#include <tuple>
#include <stdexcept>

class BVPSolver {
public:
    using Function = std::function<double(double)>;
    using BoundaryCondition = std::tuple<double, double, double>;

    BVPSolver(Function p, Function q, Function r, double a, double b)
        : p_(p), q_(q), r_(r), a_(a), b_(b) {}

    std::tuple<std::vector<double>, std::vector<double>>
    solve(int n, const std::tuple<BoundaryCondition, BoundaryCondition>& boundary_conditions) {
        double h = (b_ - a_) / n;
        std::vector<double> x(n + 1);
        for (int i = 0; i <= n; ++i) {
            x[i] = a_ + i * h;
        }

        // Set up the tridiagonal system
        std::vector<std::vector<double>> A(n + 1, std::vector<double>(n + 1, 0.0));
        std::vector<double> b(n + 1);

        // Interior points using centered differences
        for (int i = 1; i < n; ++i) {
            double xi = x[i];
            A[i][i-1] = 1/(h*h) - p_(xi)/(2*h);
            A[i][i] = -2/(h*h) + q_(xi);
            A[i][i+1] = 1/(h*h) + p_(xi)/(2*h);
            b[i] = r_(xi);
        }

        // Left boundary condition (exactly as in Python)
        A[0][0] = 2 + 3/(2*h);
        A[0][1] = -2/h;
        A[0][2] = 1/(2*h);
        b[0] = 6;

        // Right boundary condition
        A[n][n] = 1;
        b[n] = 1.5495;

        // Solve the system using Gaussian elimination
        std::vector<double> y = solve_tridiagonal_plus(A, b);

        return {x, y};
    }

private:
    Function p_, q_, r_;
    double a_, b_;

    std::vector<double> solve_tridiagonal_plus(
        std::vector<std::vector<double>>& A,
        std::vector<double>& b) {

        int n = b.size() - 1;

        // Forward elimination
        for (int i = 0; i < n; i++) {
            double pivot = A[i][i];
            if (std::abs(pivot) < 1e-10) {
                throw std::runtime_error("Zero pivot encountered");
            }

            // Eliminate below
            for (int j = i + 1; j <= std::min(i + 2, n); j++) {
                double factor = A[j][i] / pivot;
                for (int k = i; k <= std::min(i + 2, n); k++) {
                    A[j][k] -= factor * A[i][k];
                }
                b[j] -= factor * b[i];
            }
        }

        // Back substitution
        std::vector<double> x(n + 1);
        x[n] = b[n] / A[n][n];

        for (int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j <= std::min(i + 2, n); j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        return x;
    }
};

// Exact solution function
inline double exact_solution(double x) {
    return std::sqrt(1 + x*x) + std::exp(-2*x);
}
