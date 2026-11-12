import sympy as sp

# =========================================================================
# STEP 1: DEFINE VARIABLES & MONOMIAL ORDERING
# Variable precedence: x1 > y1 > z1 > x2 > y2 > z2 > x3 > y3 > z3
# =========================================================================
x1, y1, z1, x2, y2, z2, x3, y3, z3 = sp.symbols('x1 y1 z1 x2 y2 z2 x3 y3 z3')
gens = [x1, y1, z1, x2, y2, z2, x3, y3, z3]

ORDER = 'grevlex'   # 'lex', 'grlex' or 'grevlex'

# Benchmark counters (deterministic, no timing)
stats = {'div_steps': 0, 'rem_terms': 0, 'coprime_skips': 0}

def get_poly(f):
    return sp.Poly(f, *gens)

def LM_LC(f):
    """Return both Leading Monomial and Leading Coefficient in one Poly call."""
    p = sp.Poly(f, *gens)
    return p.monoms(order=ORDER)[0], p.LC(order=ORDER)

def lcm_deg(lm1, lm2):
    """Total degree of the LCM of two leading monomials (for Normal Selection Strategy)."""
    return sum(max(e1, e2) for e1, e2 in zip(lm1, lm2))

# S-Polynomial calculation for a pair (f, g)
def s_polynomial(f, g, lm_f, lc_f, lm_g, lc_g):
    lcm_mono = [max(e1, e2) for e1, e2 in zip(lm_f, lm_g)]
    m_f_exp = [l - e for l, e in zip(lcm_mono, lm_f)]
    m_g_exp = [l - e for l, e in zip(lcm_mono, lm_g)]

    m_f = sp.Mul(*[v**e for v, e in zip(gens, m_f_exp)]) / lc_f
    m_g = sp.Mul(*[v**e for v, e in zip(gens, m_g_exp)]) / lc_g
    return sp.expand(m_f * f - m_g * g)

# Multivariable polynomial long division modulo basis G
def poly_reduce(f, basis, basis_lm_lc):
    p_curr = sp.expand(f)
    rem = sp.sympify(0)
    while p_curr != 0:
        reduced = False
        lm_p, lc_p = LM_LC(p_curr)

        for g, (lm_g, lc_g) in zip(basis, basis_lm_lc):
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
# =========================================================================
R = sp.Rational
f1 = x1**2 + y1**2 + z1**2 - 21
f2 = x1**2 - 12*x1 + y1**2 + z1**2 - 9
f3 = x2**2 - 12*x2 + y2**2 + z2**2 - R(69, 5)
f4 = x2**2 - 4*x2 + y2**2 - 10*y2 + z2**2 + 11
f5 = x3**2 - 4*x3 + y3**2 - 10*y3 + z3**2 - R(271, 13)
f6 = x3**2 + y3**2 + z3**2 - R(3801, 65)
f7 = (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2 - 16
f8 = (x2 - x3)**2 + (y2 - y3)**2 + (z2 - z3)**2 - 18
f9 = (x1 - x3)**2 + (y1 - y3)**2 + (z1 - z3)**2 - 10

G = [sp.expand(f) for f in [f1, f2, f3, f4, f5, f6, f7, f8, f9]]
G_lm_lc = [LM_LC(g) for g in G]

pairs = [(i, j) for i in range(len(G)) for j in range(i + 1, len(G))]

evaluations = 0
nonzero_additions = 0
zero_reductions = 0
gen_sizes = []   # (number of terms, total degree) of each added generator

# Execute Buchberger's reduction loop with Normal Selection Strategy & 1st Criterion
while pairs:
    # Sort pairs by LCM total degree so linear reductions happen before quadratic blowup
    pairs.sort(key=lambda pair: lcm_deg(G_lm_lc[pair[0]][0], G_lm_lc[pair[1]][0]))
    i, j = pairs.pop(0)

    lm_i, lc_i = G_lm_lc[i]
    lm_j, lc_j = G_lm_lc[j]

    # Buchberger's First Criterion: Coprime leading monomials reduce to 0 automatically
    if all(min(e1, e2) == 0 for e1, e2 in zip(lm_i, lm_j)):
        stats['coprime_skips'] += 1
        zero_reductions += 1
        continue

    evaluations += 1
    spoly = s_polynomial(G[i], G[j], lm_i, lc_i, lm_j, lc_j)
    rem = poly_reduce(spoly, G, G_lm_lc)

    if rem == 0:
        zero_reductions += 1
    else:
        nonzero_additions += 1
        lm_rem, lc_rem = LM_LC(rem)
        rem_monic = sp.expand(rem / lc_rem)
        new_idx = len(G)

        rp = get_poly(rem_monic)
        gen_sizes.append((len(rp.terms()), rp.total_degree()))

        for k in range(len(G)):
            pairs.append((k, new_idx))
        G.append(rem_monic)
        G_lm_lc.append((lm_rem, sp.sympify(1)))

# Compute exact reduced Gröbner Basis using SymPy across all 9 equations and 9 variables
G_reduced = sp.groebner([f1, f2, f3, f4, f5, f6, f7, f8, f9], *gens, order=ORDER)

# =========================================================================
# STEP 3: OUTPUT RESULTS & BENCHMARK METRICS
# =========================================================================
print('=' * 65)
print('1. ALGORITHMIC STEP METRICS (' + ORDER.upper() + ' ORDERING)')
print('=' * 65)
print(f' • S-Polynomial Pairs Evaluated       : {evaluations}')
print(f' • Coprime Pairs Skipped (1st Crit.)  : {stats["coprime_skips"]}')
print(f' • Non-Zero Reductions Added          : {nonzero_additions}')
print(f' • Pairs Reduced to Zero (Total)      : {zero_reductions}')
print(f' • Total Division Steps               : {stats["div_steps"]}')
print(f' • Terms Moved to Remainder           : {stats["rem_terms"]}')
print(f' • New Generators (terms, degree)     : {gen_sizes}')

# print('\n' + '=' * 65)
# print('2. COMPUTED REDUCED GRÖBNER BASIS (' + ORDER.upper() + ' ORDERING)')
# print('=' * 65)
# for idx, poly in enumerate(G_reduced, 1):
#     print(f' g_{idx}: {poly} = 0')

# print('\n' + '=' * 65)
# print('3. 3D KINEMATIC SOLUTION EXTRACTION & FILTERING')
# print('=' * 65)

# # Solve the system directly from the reduced Gröbner basis
# solutions = sp.solve(list(G_reduced), gens, dict=True)

# # Sort solutions by z1 for stable ordering
# solutions.sort(key=lambda s: float(sp.re(sp.N(s.get(z1, 0)))))

# print(f' Total Algebraic Solutions Found: {len(solutions)}\n')

# for idx, sol in enumerate(solutions, 1):
#     if not all(v in sol for v in gens):
#         print(f'Solution {idx}: underdetermined -> {sol}\n')
#         continue

#     vals = {v: sp.N(sol[v]) for v in gens}

#     # Check if any coordinate has a non-negligible imaginary component
#     if any(abs(sp.im(x)) > 1e-9 for x in vals.values()):
#         print(f'Solution {idx}: COMPLEX (Rejected for physical 3D space)')
#         for pt_idx, (vx, vy, vz) in enumerate([(x1, y1, z1), (x2, y2, z2), (x3, y3, z3)], 1):
#             print(f'   Point {pt_idx}: ({complex(vals[vx]):.4f}, {complex(vals[vy]):.4f}, {complex(vals[vz]):.4f})')
#         print()
#         continue

# print(f'Solution {idx}: REAL PHYSICAL CONFIGURATION')
# for pt_idx, (vx, vy, vz) in enumerate([(x1, y1, z1), (x2, y2, z2), (x3, y3, z3)], 1):
#     px, py, pz = float(sp.re(vals[vx])), float(sp.re(vals[vy])), float(sp.re(vals[vz]))
#     print(f' • Point P{pt_idx} (x{pt_idx}, y{pt_idx}, z{pt_idx}) = ({px:.6f}, {py:.6f}, {pz:.6f})')
# print()