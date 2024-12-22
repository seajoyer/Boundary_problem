#include "bvp_solver.hpp"
#include <iostream>
#include <cstdlib>

int main() {
    // Coefficient functions for the differential equation
    auto p = [](double x) { return x / (1 + x*x); };
    auto q = [](double x) { return -1 / (1 + x*x); };
    auto r = [](double x) {
        return (3 - 2*x + 4*x*x) / (1 + x*x) * std::exp(-2*x);
    };

    BVPSolver solver(p, q, r, 0, 1);

    // Boundary conditions
    auto left_bc = std::make_tuple(2.0, -1.0, 6.0);   // 2u(0) - u'(0) = 6
    auto right_bc = std::make_tuple(1.0, 0.0, 1.5495); // u(1) = 1.5495
    auto boundary_conditions = std::make_tuple(left_bc, right_bc);

    // Solve with n = 20 (h = 0.05)
    int n = 20;
    auto [x, y] = solver.solve(n, boundary_conditions);

    // Create a gnuplot pipe
    FILE* gnuplot = popen("gnuplot -persist", "w");
    if (!gnuplot) {
        std::cerr << "Error opening gnuplot pipe\n";
        return 1;
    }

    // Set up gnuplot
    fprintf(gnuplot, "set title 'BVP Solution Comparison'\n");
    fprintf(gnuplot, "set xlabel 'x'\n");
    fprintf(gnuplot, "set ylabel 'y'\n");
    fprintf(gnuplot, "set grid\n");
    fprintf(gnuplot, "set key right bottom\n");  // Move legend to bottom right
    fprintf(gnuplot, "set yrange [1:2.5]\n");    // Set y-axis range to match expected solution

    // Plot exact solution first (solid blue line) and numerical solution (red dashed line with points)
    fprintf(gnuplot, "plot '-' title 'Exact solution' with lines lc rgb '#0000FF' lw 2,");
    fprintf(gnuplot, "     '-' title 'Numerical solution' with linespoints lt 2 lc rgb '#FF0000' lw 1 pt 7 ps 0.5\n");

    // Send exact solution data (more points for smooth curve)
    const int n_exact = 200;
    const double h_exact = 1.0 / (n_exact - 1);
    for (int i = 0; i < n_exact; ++i) {
        double x_exact = i * h_exact;
        fprintf(gnuplot, "%f %f\n", x_exact, exact_solution(x_exact));
    }
    fprintf(gnuplot, "e\n");

    // Send numerical solution data
    for (size_t i = 0; i < x.size(); ++i) {
        fprintf(gnuplot, "%f %f\n", x[i], y[i]);
    }
    fprintf(gnuplot, "e\n");

    pclose(gnuplot);
    return 0;
}
