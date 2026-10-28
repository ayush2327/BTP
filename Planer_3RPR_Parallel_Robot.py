import math
import sympy as sp

# =========================================================================
# STEP 1: DEFINE VARIABLES & MONOMIAL ORDERING
# Variable precedence: x > y > c > s
# =========================================================================
x, y, c, s = sp.symbols('x y c s')
gens = [x, y, c, s]          # single source of truth for the variable order

ORDER = 'lex'            # 'lex', 'grlex' or 'grevlex'

# Benchmark counters (deterministic, no timing)
stats = {'div_steps': 0, 'rem_terms': 0}

def get_poly(f):
    return sp.Poly(f, *gens)

def LM(f):
    return get_poly(f).monoms(order=ORDER)[0]

def LC(f):
    return get_poly(f).LC(order=ORDER)


SUP = str.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹')

def mono_str(exps):
    out = ''
    for v, e in zip(gens, exps):
        if e == 0:
            continue
        out += str(v) + (str(e).translate(SUP) if e > 1 else '')
    return out

def coeff_str(a):
    # a is a positive Rational; returns '' for 1 so that "1*x" prints as "x"
    if a == 1:
        return ''
    if a.q == 1:
        return str(a.p)
    return f'({a.p}/{a.q})'

def poly_str(f):
    terms = get_poly(f).terms(order=ORDER)      # leading term first
    if not terms:
        return '0'
    result = ''
    for k, (exps, coeff) in enumerate(terms):
        coeff = sp.Rational(coeff)
        neg = coeff < 0
        a = abs(coeff)
        mono = mono_str(exps)
        if mono == '':                           # constant term
            body = str(a.p) if a.q == 1 else f'{a.p}/{a.q}'
        else:
            body = coeff_str(a) + mono
        if k == 0:
            result += ('-' if neg else '') + body
        else:
            result += (' - ' if neg else ' + ') + body
    return result

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
# STEP 2: DEFINE INITIAL IDEAL GENERATORS
# =========================================================================
f1 = x**2 + y**2 + 4*x - 2*x*c - 4*c - 2*y*s + 3
f2 = x**2 + y**2 + 2*x*c - 4*x - 4*c + 2*y*s + 3
f3 = x**2 + y**2 - 2*x*s + 2*y*c - 6*y - 6*c + 9
f4 = c**2 + s**2 - 1

F = [f1, f2, f3, f4]
G = list(F)
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

# Compute exact reduced Gröbner Basis using SymPy (monic over QQ)
G_reduced = sp.groebner(F, *gens, order=ORDER, domain=sp.QQ)

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
print(f' • Final Intermediate Basis Size      : {len(G)}')

print('\n' + '=' * 65)
print('2. COMPUTED REDUCED GRÖBNER BASIS (' + ORDER.upper() + ' ORDERING)')
print('=' * 65)
for idx, poly in enumerate(G_reduced, 1):
    print(f' g{idx} = {poly_str(poly)}')

print('\n' + '=' * 65)
print('3. SOLUTION EXTRACTION')
print('=' * 65)

solutions = sp.solve(list(G_reduced), gens, dict=True)
print(f' Solutions found: {len(solutions)}\n')

valid = []
for idx, sol in enumerate(solutions, 1):
    if not all(v in sol for v in gens):
        print(f'Solution {idx}: not fully determined -> {sol}\n')
        continue

    vals = {v: sp.N(sol[v]) for v in gens}

    if any(abs(sp.im(val)) > 1e-9 for val in vals.values()):
        print(f'Solution {idx}: complex -> REJECTED')
        print('  ' + ', '.join(f'{v} = {vals[v]}' for v in gens) + '\n')
        continue

    xv, yv, cv, sv = (float(sp.re(vals[v])) for v in gens)
    theta = math.degrees(math.atan2(sv, cv))
    valid.append((xv, yv, cv, sv, theta))

valid.sort()
for idx, (xv, yv, cv, sv, theta) in enumerate(valid, 1):
    print(f'Real solution {idx}:')
    print(f'  x = {xv:.6f}, y = {yv:.6f}')
    print(f'  c = {cv:.6f}, s = {sv:.6f}  ->  theta = {theta:.2f}°\n')