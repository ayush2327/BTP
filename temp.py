
import sympy as sp
import numpy as np
from scipy.optimize import least_squares

# ==========================================
# Variables
# ==========================================

x1, y1, z1 = sp.symbols('x1 y1 z1')
x2, y2, z2 = sp.symbols('x2 y2 z2')
x3, y3, z3 = sp.symbols('x3 y3 z3')

variables = [x1, y1, z1, x2, y2, z2, x3, y3, z3]

# ==========================================
# Original polynomial generators
# ==========================================

equations = [
    x1**2 + y1**2 + z1**2 - 21,

    x1**2 - 12*x1 + y1**2 + z1**2 - 9,

    x2**2 - 12*x2 + y2**2 + z2**2 - sp.Rational(69, 5),

    x2**2 - 4*x2 + y2**2 - 10*y2 + z2**2 + 11,

    x3**2 - 4*x3 + y3**2 - 10*y3 + z3**2 - sp.Rational(271, 13),

    x3**2 + y3**2 + z3**2 - sp.Rational(3801, 65),

    (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 - 16,

    (x2-x3)**2 + (y2-y3)**2 + (z2-z3)**2 - 18,

    (x1-x3)**2 + (y1-y3)**2 + (z1-z3)**2 - 10
]

# ==========================================
# Convert equations and Jacobian to numerical functions
# ==========================================

F = sp.lambdify(variables, equations, modules="numpy")

jacobian = sp.Matrix(equations).jacobian(variables)
J = sp.lambdify(variables, jacobian, modules="numpy")


def function(v):
    return np.asarray(F(*v), dtype=float).reshape(-1)


def derivative(v):
    return np.asarray(J(*v), dtype=float)


# ==========================================
# Find distinct real solutions using multiple starts
# ==========================================

def find_real_solutions(n_starts=4000):

    rng = np.random.default_rng(42)

    solutions = []

    tolerance = 1e-8

    for i in range(n_starts):

        initial = rng.uniform(-10, 10, 9)

        result = least_squares(
            function,
            initial,
            jac=derivative,
            bounds=(-10*np.ones(9), 10*np.ones(9)),
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
            max_nfev=1500,
            x_scale="jac"
        )

        candidate = result.x

        # Check whether the candidate satisfies all equations
        error = np.linalg.norm(function(candidate), ord=np.inf)

        if error > tolerance:
            continue

        # Remove duplicate solutions
        duplicate = False

        for old in solutions:
            if np.linalg.norm(candidate - old, ord=np.inf) < 1e-6:
                duplicate = True
                break

        if not duplicate:
            solutions.append(candidate)

        if (i + 1) % 100 == 0:
            print(f"Completed {i+1}/{n_starts} initial guesses. "
                  f"Solutions found: {len(solutions)}")

    solutions.sort(key=lambda v: tuple(v))

    return solutions


# ==========================================
# Display results
# ==========================================

solutions = find_real_solutions()

np.set_printoptions(precision=10, suppress=True)

print("\n========== REAL SOLUTIONS ==========")

print("Total distinct real solutions found:", len(solutions))

for i, sol in enumerate(solutions):

    error = np.linalg.norm(function(sol), ord=np.inf)

    print(f"\nSolution {i+1}:")

    for name, value in zip(
        ["x1", "y1", "z1", "x2", "y2", "z2", "x3", "y3", "z3"],
        sol
    ):
        print(f"{name} = {value:.10f}")

    print(f"Maximum equation residual = {error:.3e}")
