"""
Prime density of Meixner-Pollaczek polynomials evaluated at primes.

MP polynomials (λ=1/4, φ=π/2):
  P_0 = 1
  P_1 = 2x
  (n+1) P_{n+1} = 2x P_n - (n - 1/2) P_{n-1}

These have RATIONAL coefficients (denominators from the (n+1) division).
We work with the "monic integer" versions: clear denominators to get
integer-coefficient polynomials, then analyse.
"""

from fractions import Fraction
from sympy import (Symbol, Poly, QQ, Rational, factor, isprime, primerange,
                   factorint, gcd as sgcd)
from math import gcd
from functools import reduce
import math

x = Symbol('x')

# ── Build MP polynomials with exact rational coefficients ────────

def mp_polys_rational(nmax):
    """Return MP polynomials P_0..P_nmax as sympy expressions with Rational coeffs."""
    polys_frac = [[Fraction(1)], [Fraction(0), Fraction(2)]]
    for n in range(1, nmax):
        pn, pm = polys_frac[n], polys_frac[n-1]
        nxt = [Fraction(0)] * (n + 2)
        for j, c in enumerate(pn):
            nxt[j+1] += 2 * c
        for j, c in enumerate(pm):
            nxt[j] -= (Fraction(n) - Fraction(1, 2)) * c
        polys_frac.append([c / (n + 1) for c in nxt])

    # Convert to sympy
    result = []
    for coeffs in polys_frac:
        expr = sum(Rational(c.numerator, c.denominator) * x**j
                   for j, c in enumerate(coeffs))
        result.append(expr)
    return result

N_POLYS = 12
mp_list = mp_polys_rational(N_POLYS)

print("MEIXNER-POLLACZEK POLYNOMIALS (λ=1/4, φ=π/2)")
print("=" * 70)
for n in range(N_POLYS):
    print(f"  P_{n}(x) = {mp_list[n]}")

# ── Clear denominators to get integer polynomials ────────────────

print(f"\n{'=' * 70}")
print("INTEGER-CLEARED VERSIONS (multiply by LCD)")
print("=" * 70)

def clear_denominators(expr):
    """Multiply by LCD of coefficients to get integer polynomial. Return (multiplier, int_poly)."""
    poly = Poly(expr, x, domain=QQ)
    coeffs = poly.all_coeffs()
    # Find LCD
    denoms = [Rational(c).q for c in coeffs]
    from math import lcm
    lcd = reduce(lcm, denoms)
    int_poly = expr * lcd
    # Simplify to remove any common factor
    int_coeffs = [int(Rational(c).p * (lcd // Rational(c).q)) for c in coeffs]
    common = reduce(gcd, [abs(c) for c in int_coeffs if c != 0])
    return lcd // common, Poly(int_poly / common, x, domain=QQ).as_expr()

int_mp = {}
for n in range(N_POLYS):
    mult, ipoly = clear_denominators(mp_list[n])
    int_mp[n] = ipoly
    print(f"  {mult}·P_{n}(x) = {ipoly}")


# ── Irreducibility ───────────────────────────────────────────────

print(f"\n{'=' * 70}")
print("IRREDUCIBILITY OVER Q")
print("=" * 70)

for n in range(N_POLYS):
    p = int_mp[n]
    if Poly(p, x).degree() < 2:
        print(f"  P_{n}: degree {Poly(p, x).degree()}, skipping")
        continue
    poly = Poly(p, x, domain=QQ)
    _, factor_pairs = poly.factor_list()
    num_factors = len(factor_pairs)
    degrees = sorted([Poly(fi, x).degree() for fi, _ in factor_pairs])
    irr = "IRREDUCIBLE" if num_factors == 1 else f"REDUCIBLE → {degrees}"
    print(f"  P_{n}: {irr}")
    if num_factors > 1:
        print(f"       = {factor(p, domain=QQ)}")


# ── Extract irreducible factors of degree ≥ 2 ───────────────────

print(f"\n{'=' * 70}")
print("IRREDUCIBLE FACTORS (degree ≥ 2)")
print("=" * 70)

irr_factors = []
seen = set()
for n in range(2, N_POLYS):
    p = int_mp[n]
    poly = Poly(p, x, domain=QQ)
    _, factor_pairs = poly.factor_list()
    for fi, mult in factor_pairs:
        deg = Poly(fi, x).degree()
        if deg >= 2:
            fi_expr = Poly(fi, x, domain=QQ).as_expr()
            fi_str = str(fi_expr)
            if fi_str not in seen:
                seen.add(fi_str)
                irr_factors.append((n, fi_expr, deg))
                print(f"  From P_{n}: {fi_expr}  (degree {deg})")


# ── Bateman-Horn constants ───────────────────────────────────────

print(f"\n{'=' * 70}")
print("BATEMAN-HORN ANALYSIS")
print("=" * 70)

sieve_primes = list(primerange(2, 200))

def count_roots_mod_p(f_expr, p):
    poly = Poly(f_expr, x, domain=QQ)
    coeffs = poly.all_coeffs()
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
        C *= (1 - w / p) / (1 - 1 / p)
    return C

# Include T_2 for comparison
from sympy import chebyshevt
T2 = chebyshevt(2, x)
T2_C = bateman_horn(T2)

print(f"\n  {'Source':>10} {'Polynomial':>35} {'deg':>4} {'C_f':>7} {'C/d':>7} {'blocked<50':>11}")
print(f"  {'-' * 80}")
blocked_T2 = sum(1 for p in primerange(2, 50) if count_roots_mod_p(T2, p) == 0)
print(f"  {'T_2':>10} {'2x²-1':>35} {2:>4} {T2_C:>7.3f} {T2_C/2:>7.4f} {blocked_T2:>11}")

for src_n, fi, deg in irr_factors:
    C = bateman_horn(fi)
    blocked = sum(1 for p in primerange(2, 50) if count_roots_mod_p(fi, p) == 0)
    short = str(fi)
    if len(short) > 35:
        short = short[:32] + "..."
    print(f"  {'MP_'+str(src_n):>10} {short:>35} {deg:>4} {C:>7.3f} {C/deg:>7.4f} {blocked:>11}")


# ── Actual prime counts at prime arguments ───────────────────────

print(f"\n{'=' * 70}")
print("ACTUAL PRIME COUNTS: |f(p)| for primes p < 500")
print("=" * 70)

eval_primes = list(primerange(2, 500))

# First: raw MP polynomials
print("\nRaw MP polynomials P_n(p):")
for n in range(2, min(N_POLYS, 9)):
    P = mp_list[n]
    count = 0
    for p in eval_primes:
        val = P.subs(x, p)
        # Check if it's an integer
        if val == int(val):
            val = abs(int(val))
            if val > 1 and isprime(val):
                count += 1
    print(f"  P_{n}(p): {count}/{len(eval_primes)} prime")

# Irreducible factors
print("\nIrreducible factors:")
for src_n, fi, deg in irr_factors:
    C = bateman_horn(fi)
    count = 0
    first_primes = []
    for p in eval_primes:
        val = abs(int(fi.subs(x, p)))
        if val > 1 and isprime(val):
            count += 1
            if len(first_primes) < 5:
                first_primes.append((p, val))
    density = count / len(eval_primes)
    print(f"  From MP_{src_n} (deg {deg}): {count}/{len(eval_primes)} = {density:.3f}  "
          f"(C_f/deg = {C/deg:.4f})")
    if first_primes:
        print(f"    First hits: {first_primes}")


# ── Shield comparison ────────────────────────────────────────────

print(f"\n{'=' * 70}")
print("SHIELD: which small primes are blocked?")
print("=" * 70)

all_polys = [("T_2", T2, 2)] + [(f"MP_{n}", fi, deg) for n, fi, deg in irr_factors if deg <= 4]

print(f"\n  {'q':>4}", end="")
for name, _, _ in all_polys:
    print(f"{name:>10}", end="")
print()

for q in primerange(2, 50):
    print(f"  {q:>4}", end="")
    for name, fi, _ in all_polys:
        w = count_roots_mod_p(fi, q)
        if w == 0:
            print(f"{'·':>10}", end="")
        else:
            print(f"{w:>10}", end="")
    print()


# ── MP_2 vs T_2: direct comparison ──────────────────────────────

print(f"\n{'=' * 70}")
print("HEAD-TO-HEAD: T_2(x) = 2x²-1 vs MP_2 irreducible factor")
print("=" * 70)

if irr_factors:
    mp2_src, mp2_poly, mp2_deg = irr_factors[0]  # first degree-2 factor
    if mp2_deg == 2:
        print(f"\n  T_2(x) = {T2}")
        print(f"  MP factor = {mp2_poly}")
        print(f"\n  {'p':>4} {'T_2(p)':>10} {'prime?':>7} {'MP(p)':>10} {'prime?':>7}")
        print(f"  {'-' * 42}")
        for p in list(primerange(2, 50)):
            t2_val = int(T2.subs(x, p))
            mp_val = abs(int(mp2_poly.subs(x, p)))
            t2_prime = "✓" if isprime(abs(t2_val)) else ""
            mp_prime = "✓" if isprime(mp_val) else ""
            print(f"  {p:>4} {t2_val:>10} {t2_prime:>7} {mp_val:>10} {mp_prime:>7}")
