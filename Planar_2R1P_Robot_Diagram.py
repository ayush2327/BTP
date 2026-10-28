import math
import sympy as sp

# =========================================================================
# STEP 1: DEFINE VARIABLES & MONOMIAL ORDERING
# Variable precedence: c1 > s1 > c2 > s2 > p
# =========================================================================
c1, s1, c2, s2, p = sp.symbols('c1 s1 c2 s2 p')
gens = [c1, s1, c2, s2, p]

def get_poly(f):
    return sp.Poly(f, c1, s1, c2, s2, p)

ORDER = 'grevlex'   # 'lex', 'grlex' or 'grevlex'

# Physical limits of the prismatic joint (problem-specific)
P_MIN, P_MAX = 0.0, 1.0

# Benchmark counters (deterministic, no timing)
stats = {'div_steps': 0, 'rem_terms': 0}

def LM(f):
    return get_poly(f).monoms(order=ORDER)[0]

def LC(f):
    return get_poly(f).LC(order=ORDER)

# S-Polynomial calculation for a pair (f, g)
def s_polynomial(f, g):
    lm_f, lm_g = LM(f), LM(g)
    lc_f, lc_g = LC(f), LC(g)

    # Least Common Multiple (LCM) of leading monomials
    lcm_mono = [max(e1, e2) for e1, e2 in zip(lm_f, lm_g)]
    m_f_exp = [l - e for l, e in zip(lcm_mono, lm_f)]
    m_g_exp = [l - e for l, e in zip(lcm_mono, lm_g)]

    m_f = sp.Mul(*[v**e for v, e in zip(gens, m_f_exp)]) / lc_f
    m_g = sp.Mul(*[v**e for v, e in zip(gens, m_g_exp)]) / lc_g
    return sp.expand(m_f * f - m_g * g)

# Multivariable polynomial long division modulo basis G
def poly_reduce(f, basis):
    p_curr = sp.expand(f)
    rem = sp.sympify(0)
    while p_curr != 0:
        reduced = False
        lm_p, lc_p = LM(p_curr), LC(p_curr)

        for g in basis:
            lm_g, lc_g = LM(g), LC(g)

            # Check if lm_g divides lm_p
            if all(e_p >= e_g for e_p, e_g in zip(lm_p, lm_g)):
                stats['div_steps'] += 1
                m_exp = [e_p - e_g for e_p, e_g in zip(lm_p, lm_g)]
                m = (lc_p / lc_g) * sp.Mul(*[v**e for v, e in zip(gens, m_exp)])
                p_curr = sp.expand(p_curr - m * g)
                reduced = True
                break

        if not reduced:
            stats['rem_terms'] += 1
            lt_term = lc_p * sp.Mul(*[v**e for v, e in zip(gens, lm_p)])
            rem += lt_term
            p_curr = sp.expand(p_curr - lt_term)
    return rem

# =========================================================================
# STEP 2: DEFINE INITIAL KINEMATIC IDEAL GENERATORS
# (Target: a = 3/2, b = 39/20, theta = 0 deg)
# =========================================================================
f1 = c1**2 + s1**2 - 1
f2 = c2**2 + s2**2 - 1
f3 = c1 * c2 - s1 * s2 - 1
f4 = c1 * s2 + c2 * s1
f5 = 2 * c1 + p - sp.Rational(1, 2)
f6 = 2 * s1 - sp.Rational(39, 20)

G = [f1, f2, f3, f4, f5, f6]
pairs = [(i, j) for i in range(len(G)) for j in range(i + 1, len(G))]

evaluations = 0
nonzero_additions = 0
zero_reductions = 0
gen_sizes = []   # (number of terms, total degree) of each added generator

# Execute Buchberger's reduction loop
while pairs:
    i, j = pairs.pop(0)
    evaluations += 1
    spoly = s_polynomial(G[i], G[j])
    rem = poly_reduce(spoly, G)

    if rem == 0:
        zero_reductions += 1
    else:
        nonzero_additions += 1
        lc = LC(rem)
        rem_monic = sp.expand(rem / lc)
        new_idx = len(G)

        rp = get_poly(rem_monic)
        gen_sizes.append((len(rp.terms()), rp.total_degree()))

        for k in range(len(G)):
            pairs.append((k, new_idx))
        G.append(rem_monic)

# Compute exact reduced Gröbner Basis using SymPy
G_reduced = sp.groebner([f1, f2, f3, f4, f5, f6], c1, s1, c2, s2, p, order=ORDER)

# =========================================================================
# STEP 3: OUTPUT RESULTS & BENCHMARK METRICS
# =========================================================================
print('=' * 65)
print('1. ALGORITHMIC STEP METRICS (' + ORDER.upper() + ' ORDERING)')
print('=' * 65)
print(f' • Total S-Polynomial Pairs Evaluated : {evaluations}')
print(f' • Non-Zero Reductions Added          : {nonzero_additions}')
print(f' • Pairs Reduced to Zero              : {zero_reductions}')
print(f' • Total Division Steps               : {stats["div_steps"]}')
print(f' • Terms Moved to Remainder           : {stats["rem_terms"]}')
print(f' • New Generators (terms, degree)     : {gen_sizes}')

print('\n' + '=' * 65)
print('2. COMPUTED REDUCED GRÖBNER BASIS (' + ORDER.upper() + ' ORDERING)')
print('=' * 65)
for idx, poly in enumerate(G_reduced, 1):
    print(f' g_{idx}: {poly} = 0')

print('\n' + '=' * 65)
print('3. KINEMATIC SOLUTION EXTRACTION & PHYSICAL FILTERING')
print('=' * 65)

# Solve the system directly from the Gröbner basis (no hardcoded values)
solutions = sp.solve(list(G_reduced), gens, dict=True)

# Sort by p so the output order is stable across monomial orderings
solutions.sort(key=lambda s: float(sp.re(sp.N(s.get(p, 0)))))

print(f' Solutions found: {len(solutions)}\n')

for idx, sol in enumerate(solutions, 1):
    # Make sure every variable was determined (zero-dimensional system)
    if not all(v in sol for v in gens):
        print(f'Solution {idx}: not fully determined -> {sol}\n')
        continue

    vals = {v: sp.N(sol[v]) for v in gens}

    # Reject complex solutions
    if any(abs(sp.im(x)) > 1e-9 for x in vals.values()):
        print(f'Solution {idx}: complex -> REJECTED\n')
        continue

    c1_val, s1_val = float(sp.re(vals[c1])), float(sp.re(vals[s1]))
    c2_val, s2_val = float(sp.re(vals[c2])), float(sp.re(vals[s2]))
    p_val = float(sp.re(vals[p]))

    print(f'Solution {idx}: p = {p_val:.6f}')

    # Physical filter on the prismatic joint
    if P_MIN <= p_val <= P_MAX:
        theta1_deg = math.degrees(math.atan2(s1_val, c1_val))
        theta2_deg = math.degrees(math.atan2(s2_val, c2_val))

        print(f' -> VALID: Satisfies physical limits ({P_MIN} <= p <= {P_MAX})')
        print(f' c1 = {c1_val:.4f}, s1 = {s1_val:.4f}, c2 = {c2_val:.4f}, s2 = {s2_val:.4f}')
        print(f' • Revolute Joint 1 (theta1) = {theta1_deg:.2f}°')
        print(f' • Revolute Joint 2 (theta2) = {theta2_deg:.2f}°')
        print(f' • Prismatic Extension (p) = {p_val:.4f} units\n')
    else:
        print(f' -> REJECTED: Outside physical slider range [{P_MIN} <= p <= {P_MAX}]\n')