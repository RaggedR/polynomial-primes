"""
Bateman-Horn analysis for Hermite polynomials evaluated at primes.

H_0(x) = 1
H_1(x) = 2x
H_2(x) = 4x² - 2
H_3(x) = 8x³ - 12x
H_4(x) = 16x⁴ - 48x² + 12
H_5(x) = 32x⁵ - 160x³ + 120x
H_6(x) = 64x⁶ - 480x⁴ + 720x² - 120

But first: are they even irreducible?
"""

from sympy import hermite, Symbol, Poly, QQ, factor, isprime, primerange, factorint
import math

x = Symbol('x')

# ── Irreducibility check ─────────────────────────────────────────

print("=" * 70)
print("HERMITE POLYNOMIALS: IRREDUCIBILITY OVER Q")
print("=" * 70)

for n in range(1, 13):
    H = hermite(n, x)
    poly = Poly(H, x, domain=QQ)
    _, factor_pairs = poly.factor_list()
    num_factors = len(factor_pairs)
    degrees = sorted([Poly(fi, x).degree() for fi, _ in factor_pairs])
    irr = "IRREDUCIBLE" if num_factors == 1 else f"REDUCIBLE → {degrees}"
    print(f"  H_{n:>2}(x) = {H}")
    print(f"         {irr}")
    if num_factors > 1:
        print(f"         = {factor(H, domain=QQ)}")
    print()


# ── Factor out common structure ──────────────────────────────────

print("=" * 70)
print("STRUCTURAL FACTORS")
print("=" * 70)
print("""
For odd n:  H_n(x) = x · Q(x²)     → always has factor x → REDUCIBLE
For even n: H_n(x) = R(x²) + const → polynomial in x², might be irreducible
""")

# The irreducible pieces are what matter for Bunyakovsky
print("Extracting irreducible factors of H_n:")
for n in range(1, 13):
    H = hermite(n, x)
    poly = Poly(H, x, domain=QQ)
    content, factor_pairs = poly.factor_list()
    for fi, mult in factor_pairs:
        deg = Poly(fi, x).degree()
        if deg >= 2:
            print(f"  H_{n:>2} irreducible factor (deg {deg}): {fi}")


# ── Bateman-Horn for each irreducible factor ─────────────────────

print(f"\n{'=' * 70}")
print("BATEMAN-HORN CONSTANTS FOR IRREDUCIBLE FACTORS")
print("=" * 70)

sieve_primes = list(primerange(2, 200))

def count_roots_mod_p(f_poly, p):
    coeffs = Poly(f_poly, x, domain=QQ).all_coeffs()
    count = 0
    for x_val in range(p):
        val = 0
        for c in coeffs:
            val = (val * x_val + int(c)) % p
        if val == 0:
            count += 1
    return count

def bateman_horn(f_poly):
    C = 1.0
    for p in sieve_primes:
        w = count_roots_mod_p(f_poly, p)
        factor = (1 - w / p) / (1 - 1 / p)
        C *= factor
    return C

# Collect all irreducible factors of degree ≥ 2
irr_factors = []
for n in range(1, 13):
    H = hermite(n, x)
    poly = Poly(H, x, domain=QQ)
    content, factor_pairs = poly.factor_list()
    for fi, mult in factor_pairs:
        deg = Poly(fi, x).degree()
        if deg >= 2:
            # Check content = gcd of coefficients
            fi_poly = Poly(fi, x, domain=QQ)
            # Check fixed divisor
            vals = [abs(int(fi_poly.eval(v))) for v in range(deg + 1)]
            from math import gcd
            from functools import reduce
            fixed_div = reduce(gcd, [v for v in vals if v != 0])
            irr_factors.append((n, fi, deg, fixed_div))

print(f"\n  {'Source':>8} {'Degree':>7} {'Fixed div':>10} {'C_f':>8} {'C_f/deg':>8} {'Polynomial'}")
print(f"  {'-' * 75}")

# Also include T_2 for comparison
T2 = 2*x**2 - 1
C_T2 = bateman_horn(T2)
print(f"  {'T_2':>8} {2:>7} {1:>10} {C_T2:>8.3f} {C_T2/2:>8.4f} {T2}")

for n, fi, deg, fd in irr_factors:
    C = bateman_horn(fi)
    print(f"  {'H_'+str(n):>8} {deg:>7} {fd:>10} {C:>8.3f} {C/deg:>8.4f} {fi}")


# ── Actual prime counts ──────────────────────────────────────────

print(f"\n{'=' * 70}")
print("ACTUAL PRIME COUNTS: |H_n(p)| for primes p < N")
print("=" * 70)

eval_primes = list(primerange(2, 500))

for n in range(1, 9):
    H = hermite(n, x)
    count = 0
    total = 0
    for p in eval_primes:
        val = abs(int(H.subs(x, p)))
        total += 1
        if val > 1 and isprime(val):
            count += 1
    density = count / total if total else 0
    print(f"  H_{n}(p) prime for {count}/{total} primes p < 500  "
          f"(density {density:.3f})")


# ── Relationship between H_n and T_n ────────────────────────────

print(f"\n{'=' * 70}")
print("H_n vs T_n: ARE THEY RELATED?")
print("=" * 70)

print("""
H_2(x) = 4x² - 2 = 2(2x² - 1) = 2·T_2(x)
So |H_2(p)| is prime iff T_2(p) = (H_2(p))/2 is prime...
but H_2(p) = 2·T_2(p) is always even, so H_2(p) itself is NEVER prime (for p≥2).

What about the irreducible factor of H_2 after removing the 2?
That's just T_2(x) = 2x² - 1 again.
""")

# Check: is the irreducible core of H_n related to T_n?
print("Irreducible cores of H_n vs Chebyshev:")
for n in range(2, 9):
    H = hermite(n, x)
    poly = Poly(H, x, domain=QQ)
    content, factor_pairs = poly.factor_list()
    cores = [(fi, Poly(fi,x).degree()) for fi, _ in factor_pairs if Poly(fi,x).degree() >= 2]
    print(f"  H_{n}: content={content}, cores: {[(str(fi), d) for fi, d in cores]}")


# ── The right question: irreducible factors evaluated at primes ──

print(f"\n{'=' * 70}")
print("PRIME DENSITY OF IRREDUCIBLE FACTORS AT PRIME ARGUMENTS")
print("=" * 70)

results = []
seen = set()

for n in range(2, 10):
    H = hermite(n, x)
    poly = Poly(H, x, domain=QQ)
    content, factor_pairs = poly.factor_list()
    for fi, mult in factor_pairs:
        deg = Poly(fi, x).degree()
        if deg < 2:
            continue
        fi_str = str(fi)
        if fi_str in seen:
            continue
        seen.add(fi_str)

        C = bateman_horn(fi)
        count = 0
        total = len(eval_primes)
        for p in eval_primes:
            val = abs(int(Poly(fi, x, domain=QQ).eval(p)))
            if val > 1 and isprime(val):
                count += 1
        density = count / total

        # Blocked primes
        blocked = [p for p in primerange(2, 50) if count_roots_mod_p(fi, p) == 0]

        results.append((f"H_{n}", fi, deg, C, C/deg, density, blocked))
        print(f"\n  From H_{n}: f(x) = {fi}  (degree {deg})")
        print(f"    C_f = {C:.3f},  C_f/deg = {C/deg:.4f}")
        print(f"    Actual density at primes p<500: {count}/{total} = {density:.3f}")
        print(f"    Blocked primes < 50: {blocked}")

# Comparison summary
print(f"\n{'=' * 70}")
print("COMPARISON: CHEBYSHEV vs HERMITE IRREDUCIBLE FACTORS")
print("=" * 70)
print(f"\n  {'Source':>10} {'Poly':>25} {'deg':>4} {'C_f':>7} {'C/d':>7} {'density':>8} {'blocked':>8}")
print(f"  {'-' * 75}")
print(f"  {'T_2':>10} {'2x²-1':>25} {2:>4} {C_T2:>7.3f} {C_T2/2:>7.4f} {'':>8} {9:>8}")
for src, fi, deg, C, Cd, dens, bl in results:
    print(f"  {src:>10} {str(fi):>25} {deg:>4} {C:>7.3f} {Cd:>7.4f} {dens:>8.3f} {len(bl):>8}")
