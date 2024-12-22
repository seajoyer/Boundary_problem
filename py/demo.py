import numpy as np
import matplotlib.pyplot as plt
from bvp_solver import BVPSolver, exact_solution

def main():
    # Coefficient functions for the differential equation
    def p(x):
        return x / (1 + x**2)

    def q(x):
        return -1 / (1 + x**2)

    def r(x):
        return (3 - 2*x + 4*x**2) / (1 + x**2) * np.exp(-2*x)

    solver = BVPSolver(p, q, r, 0, 1)

    # Boundary conditions
    left_bc = (2, -1, 6)       # 2u(0) - u'(0) = 6
    right_bc = (1, 0, 1.5495)  # u(1) = 1.5495
    boundary_conditions = (left_bc, right_bc)

    n = 20  # h = 0.05 = 1/20
    x, y = solver.solve(n, boundary_conditions)

    x_exact = np.linspace(0, 1, 200)
    y_exact = exact_solution(x_exact)

    plt.figure(figsize=(10, 6))
    plt.plot(x_exact, y_exact, 'b-', label='Exact solution')
    plt.plot(x, y, 'ro--', label='Numerical solution', markersize=4)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('BVP Solution Comparison')
    plt.grid(True)
    plt.legend()

    plt.show()

if __name__ == "__main__":
    main()
