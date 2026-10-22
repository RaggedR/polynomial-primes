"""
Survey of orthogonal polynomial families for prime-friendliness.

For each family, extract irreducible factors, compute Bateman-Horn C_f,
and measure actual prime density at prime arguments.

Families:
  1. Chebyshev T_n (done — baseline)
  2. Chebyshev U_n (second kind)
  3. Legendre P_n
  4. Laguerre L_n
  5. Gegenbauer/ultraspherical C_n^λ
  6. Dickson D_n(x, a) for a ≠ 1
  7. Bernoulli B_n(x)
  8. Euler E_n(x)
  9. Cyclotomic Φ_n(x)
  10. Fibonacci polynomials
  11. Krawtchouk (discrete)
  12. Charlier (discrete)
"""

from sympy import (
    Symbol, Poly, QQ, Rational, factor, isprime, primerange,
    chebyshevt, chebyshevu, legendre_poly, laguerre, bernoulli,
    euler, gegenbauer, cyclotomic_poly, fibonacci,
    sqrt, binomial, factorial, symbols
)
from functools import reduce
from math import gcd
import math

x = Symbol('x')
n_sym = Symbol('n')

sieve_primes = list(primerange(2, 200))
eval_primes = list(primerange(2, 2000))

def count_roots_mod_p(f_expr, p):
    try:
        poly = Poly(f_expr, x, domain=QQ)
        coeffs = poly.all_coeffs()
    except Exception:
        return -1
    count = 0
    for x_val in range(p):
        val = 0
        for c in coeffs:
            val = (val * x_val + int(c)) % p
        if val == 0:
            count += 1
    return count

def bateman_horn(f_expr):
    C = 1.0
    for p in sieve_primes:
        w = count_roots_mod_p(f_expr, p)
        if w < 0:
            return -1
        f = (1 - w / p) / (1 - 1 / p)
        C *= f
        if C < 1e-10:
            return 0.0
    return C

def analyse_poly(name, f_expr):
    """Analyse a single polynomial: irreducibility, C_f, prime density."""
    try:
        poly = Poly(f_expr, x, domain=QQ)
    except Exception:
        return None

    deg = poly.degree()
    if deg < 2:
        return None

    # Factor over Q
    _, factor_pairs = poly.factor_list()

    results = []
    seen = set()
    for fi, mult in factor_pairs:
        d = Poly(fi, x).degree()
        if d < 2:
            continue
        fi_expr = Poly(fi, x, domain=QQ).as_expr()
        fi_str = str(fi_expr)
        if fi_str in seen:
            continue
        seen.add(fi_str)

        C = bateman_horn(fi_expr)
        if C < 0:
            continue

        # Prime density at prime arguments
        count = 0
        for p in eval_primes:
            try:
                val = abs(int(fi_expr.subs(x, p)))
                if val > 1 and isprime(val):
                    count += 1
            except Exception:
                pass
        density = count / len(eval_primes)
        blocked = sum(1 for p in primerange(2, 50) if count_roots_mod_p(fi_expr, p) == 0)

        results.append({
            'source': name,
            'poly': fi_expr,
            'deg': d,
            'C_f': C,
            'C_d': C / d if d > 0 else 0,
            'density': density,
            'blocked': blocked,
            'n_irr': len(factor_pairs)
        })

    return results


# ══════════════════════════════════════════════════════════════════
# BUILD ALL POLYNOMIAL FAMILIES
# ══════════════════════════════════════════════════════════════════

all_results = []

# ── 1. Chebyshev T_n (baseline) ──────────────────────────────────
print("Chebyshev T_n...")
for n in range(2, 13):
    r = analyse_poly(f"T_{n}", chebyshevt(n, x))
    if r: all_results.extend(r)

# ── 2. Chebyshev U_n (second kind) ───────────────────────────────
print("Chebyshev U_n...")
for n in range(2, 13):
    r = analyse_poly(f"U_{n}", chebyshevu(n, x))
    if r: all_results.extend(r)

# ── 3. Legendre P_n ──────────────────────────────────────────────
# Legendre has rational coefficients; clear denominators
print("Legendre P_n...")
for n in range(2, 13):
    P = legendre_poly(n, x)
    # Clear denominators
    poly = Poly(P, x, domain=QQ)
    coeffs = poly.all_coeffs()
    denoms = [Rational(c).q for c in coeffs]
    from math import lcm
    lcd = reduce(lcm, denoms)
    int_poly = P * lcd
    # Remove common factor
    int_coeffs = [abs(int(c * lcd)) for c in coeffs]
    common = reduce(gcd, [c for c in int_coeffs if c != 0])
    int_poly = int_poly / common
    r = analyse_poly(f"Leg_{n}", int_poly)
    if r: all_results.extend(r)

# ── 4. Laguerre L_n ──────────────────────────────────────────────
print("Laguerre L_n...")
for n in range(2, 10):
    L = laguerre(n, x)
    poly = Poly(L, x, domain=QQ)
    coeffs = poly.all_coeffs()
    denoms = [Rational(c).q for c in coeffs]
    lcd = reduce(lcm, denoms)
    int_poly = L * lcd
    int_coeffs = [abs(int(c * lcd)) for c in coeffs]
    common = reduce(gcd, [c for c in int_coeffs if c != 0])
    int_poly = int_poly / common
    r = analyse_poly(f"Lag_{n}", int_poly)
    if r: all_results.extend(r)

# ── 5. Gegenbauer C_n^λ for various λ ────────────────────────────
print("Gegenbauer C_n^λ...")
for lam_num, lam_den in [(1, 2), (1, 3), (2, 3), (1, 4), (3, 4), (1, 1), (2, 1), (3, 2)]:
    lam = Rational(lam_num, lam_den)
    for n in range(2, 8):
        try:
            G = gegenbauer(n, lam, x)
            poly = Poly(G, x, domain=QQ)
            coeffs = poly.all_coeffs()
            denoms = [Rational(c).q for c in coeffs]
            lcd = reduce(lcm, denoms)
            int_poly = G * lcd
            int_coeffs = [abs(int(c * lcd)) for c in coeffs]
            common = reduce(gcd, [c for c in int_coeffs if c != 0])
            int_poly = int_poly / common
            r = analyse_poly(f"Geg({lam})_{n}", int_poly)
            if r: all_results.extend(r)
        except Exception:
            pass

# ── 6. Dickson D_n(x, a) for a = 2, 3, -1 ───────────────────────
print("Dickson D_n(x, a)...")
# D_n(x, a): D_0=2, D_1=x, D_n = x·D_{n-1} - a·D_{n-2}
def dickson(n, a_val):
    if n == 0: return 2 + 0*x
    if n == 1: return x
    d0, d1 = 2 + 0*x, x
    for _ in range(2, n + 1):
        d0, d1 = d1, x * d1 - a_val * d0
    return d1.expand()

for a_val in [2, 3, -1, 5, 7]:
    for n in range(2, 10):
        D = dickson(n, a_val)
        r = analyse_poly(f"D_{n}(a={a_val})", D)
        if r: all_results.extend(r)

# ── 7. Cyclotomic Φ_n(x) ────────────────────────────────────────
print("Cyclotomic Φ_n...")
for n in range(2, 31):
    Phi = cyclotomic_poly(n, x)
    r = analyse_poly(f"Φ_{n}", Phi)
    if r: all_results.extend(r)

# ── 8. Bernoulli B_n(x) ─────────────────────────────────────────
print("Bernoulli B_n...")
for n in range(2, 10):
    B = bernoulli(n, x)
    poly = Poly(B, x, domain=QQ)
    coeffs = poly.all_coeffs()
    denoms = [Rational(c).q for c in coeffs]
    lcd = reduce(lcm, denoms)
    int_poly = B * lcd
    int_coeffs = [abs(int(c * lcd)) for c in coeffs]
    common = reduce(gcd, [c for c in int_coeffs if c != 0])
    if common > 0:
        int_poly = int_poly / common
    r = analyse_poly(f"Bern_{n}", int_poly)
    if r: all_results.extend(r)

# ── 9. Euler E_n(x) ──────────────────────────────────────────────
print("Euler E_n...")
for n in range(2, 10):
    try:
        E = euler(n, x)
        poly = Poly(E, x, domain=QQ)
        coeffs = poly.all_coeffs()
        denoms = [Rational(c).q for c in coeffs]
        lcd = reduce(lcm, denoms)
        int_poly = E * lcd
        int_coeffs = [abs(int(c * lcd)) for c in coeffs]
        common = reduce(gcd, [c for c in int_coeffs if c != 0])
        if common > 0:
            int_poly = int_poly / common
        r = analyse_poly(f"Euler_{n}", int_poly)
        if r: all_results.extend(r)
    except Exception:
        pass

# ── 10. Fibonacci polynomials ────────────────────────────────────
print("Fibonacci polynomials...")
# F_1 = 1, F_2 = x, F_{n+1} = x·F_n + F_{n-1}
def fib_poly(n):
    if n == 1: return 1 + 0*x
    if n == 2: return x
    f0, f1 = 1 + 0*x, x
    for _ in range(3, n + 1):
        f0, f1 = f1, (x * f1 + f0).expand()
    return f1

for n in range(3, 13):
    F = fib_poly(n)
    r = analyse_poly(f"Fib_{n}", F)
    if r: all_results.extend(r)


# ══════════════════════════════════════════════════════════════════
# DEDUPLICATE AND RANK
# ══════════════════════════════════════════════════════════════════

print(f"\n{'=' * 80}")
print("COMPLETE RANKING: ALL IRREDUCIBLE FACTORS BY C_f/degree")
print("=" * 80)

# Deduplicate by polynomial string
unique = {}
for r in all_results:
    key = str(r['poly'])
    if key not in unique or r['C_d'] > unique[key]['C_d']:
        unique[key] = r

ranked = sorted(unique.values(), key=lambda r: -r['C_d'])

print(f"\n  {'Rank':>4} {'Source':>15} {'deg':>4} {'C_f':>7} {'C_f/deg':>8} "
      f"{'density':>8} {'blocked':>8}  Polynomial")
print(f"  {'-' * 95}")

for i, r in enumerate(ranked[:40]):
    poly_str = str(r['poly'])
    if len(poly_str) > 35:
        poly_str = poly_str[:32] + "..."
    marker = " ★" if r['C_d'] > 1.5 else ""
    print(f"  {i+1:>4} {r['source']:>15} {r['deg']:>4} {r['C_f']:>7.3f} {r['C_d']:>8.4f} "
          f"{r['density']:>8.3f} {r['blocked']:>8}  {poly_str}{marker}")

# ── Highlight the top degree-2 polynomials ───────────────────────

print(f"\n{'=' * 80}")
print("TOP 15 DEGREE-2 POLYNOMIALS (the sweet spot)")
print("=" * 80)

deg2 = [r for r in ranked if r['deg'] == 2]
print(f"\n  {'Rank':>4} {'Source':>15} {'C_f':>7} {'C_f/deg':>8} "
      f"{'density':>8} {'blocked<50':>11}  Polynomial")
print(f"  {'-' * 80}")

for i, r in enumerate(deg2[:15]):
    print(f"  {i+1:>4} {r['source']:>15} {r['C_f']:>7.3f} {r['C_d']:>8.4f} "
          f"{r['density']:>8.3f} {r['blocked']:>11}  {r['poly']}")
