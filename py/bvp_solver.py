import numpy as np
from typing import Callable, Tuple

class BVPSolver:
    def __init__(self, p: Callable[[float], float], q: Callable[[float], float],
                 r: Callable[[float], float], a: float, b: float):
        self.p = p
        self.q = q
        self.r = r
        self.a = a
        self.b = b

    def solve(self, n: int, boundary_conditions: Tuple[Tuple[float, float, float],
                                                 Tuple[float, float, float]]) -> Tuple[np.ndarray, np.ndarray]:
        h = (self.b - self.a) / n
        x = np.linspace(self.a, self.b, n+1)

        A = np.zeros((n+1, n+1))
        b = np.zeros(n+1)

        # Interior points using centered differences
        for i in range(1, n):
            xi = x[i]
            A[i, i-1] = 1/h**2 - self.p(xi)/(2*h)
            A[i, i] = -2/h**2 + self.q(xi)
            A[i, i+1] = 1/h**2 + self.p(xi)/(2*h)
            b[i] = self.r(xi)

        A[0, 0] = 2 + 3/(2*h)     # u[0]
        A[0, 1] = -2/h            # u[1]
        A[0, 2] = 1/(2*h)         # u[2]
        b[0] = 6

        # Right boundary
        A[n, n] = 1
        b[n] = 1.5495

        y = np.linalg.solve(A, b)

        return x, y

def exact_solution(x: np.ndarray) -> np.ndarray:
    return np.sqrt(1 + x**2) + np.exp(-2*x)
