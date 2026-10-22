"""
Can we beat Euler's p²+p+41 (C_f/deg = 3.31) with a degree-5 polynomial?

We need C_f/deg > 3.31, so C_f > 16.55 for a quintic.
That means we need an ENORMOUS shield — blocking almost all small primes.

Strategy:
  1. Search for quintics with many blocked small primes
  2. The key insight: a polynomial f blocks prime q iff f(x) ≡ 0 mod q
     has no solutions. For a degree-5 poly mod q, there are q^5 possible
     coefficient tuples — we need to find ones where ω(q) = 0 for many q.
  3. We can construct candidates by taking products of known shield-heavy
     polynomials, or by direct search optimising C_f.
"""

from sympy import isprime, primerange, Poly, QQ, Symbol, factor, Rational
from functools import reduce
from math import gcd, log
import itertools
import time

x = Symbol('x')
sieve_primes = list(primerange(2, 500))
eval_primes = list(primerange(2, 10000))

def count_roots_mod_p(coeffs, p):
    """coeffs = [a5, a4, a3, a2, a1, a0] for a5*x^5 + ... + a0."""
    count = 0
    for xv in range(p):
        val = 0
        for c in coeffs:
            val = (val * xv + c) % p
        if val == 0:
            count += 1
    return count

def bateman_horn_fast(coeffs, primes=None):
    if primes is None:
        primes = sieve_primes
    C = 1.0
    for p in primes:
        w = count_roots_mod_p(coeffs, p)
        C *= (1 - w/p) / (1 - 1/p)
        if C < 1e-15:
            return 0.0
    return C

def eval_poly(coeffs, val):
    result = 0
    for c in coeffs:
        result = result * val + c
    return result

def prime_density(coeffs):
    count = 0
    for p in eval_primes:
        val = abs(eval_poly(coeffs, p))
        if val > 1 and isprime(val):
            count += 1
    return count, count / len(eval_primes)

# ── Euler's baseline ─────────────────────────────────────────────

euler_C = bateman_horn_fast([1, 1, 41])
euler_Cd = euler_C / 2
print(f"EULER BASELINE: p²+p+41")
print(f"  C_f = {euler_C:.3f}, C_f/deg = {euler_Cd:.4f}")
print(f"  TARGET for quintic: C_f > {euler_Cd * 5:.2f}")

# ── Strategy 1: Cyclotomic-inspired quintics ─────────────────────
# Φ_11(x) = x^10 + x^9 + ... + x + 1 has factors over Q.
# Φ_11(x) = (x^5 + x^4 - 4x^3 + ... ) etc? No, Φ_11 is irreducible of degree 10.
# But Φ_p for prime p gives degree p-1, so Φ_6 = x²-x+1 (degree 2)...
# No degree-5 cyclotomic polynomial exists (φ(n) = 5 has no solution).

# ── Strategy 2: Brute force search over small coefficients ───────

print(f"\n{'='*70}")
print("SEARCH 1: QUINTICS x^5 + ax^4 + bx^3 + cx^2 + dx + e")
print(f"{'='*70}")

# For speed, first screen on small primes only
screen_primes = list(primerange(2, 100))
best = []

t0 = time.time()

# Search x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e
# with small |coefficients|
for a in range(-5, 6):
    for b in range(-20, 21, 1):
        for c in range(-40, 41, 2):
            for d in range(-60, 61, 3):
                for e in [1, -1, 2, -2, 3, -3, 5, -5, 7, -7,
                          11, -11, 13, -13, 17, -17, 41, -41]:
                    coeffs = [1, a, b, c, d, e]
                    # Quick screen: how many of first 15 primes are blocked?
                    blocked = sum(1 for p in screen_primes[:15]
                                  if count_roots_mod_p(coeffs, p) == 0)
                    if blocked < 10:
                        continue
                    # Deeper check
                    C = bateman_horn_fast(coeffs, screen_primes)
                    if C / 5 > euler_Cd * 0.8:  # within 80% of target
                        # Full computation
                        C_full = bateman_horn_fast(coeffs)
                        if C_full / 5 > euler_Cd * 0.7:
                            best.append((C_full, C_full/5, coeffs))

elapsed = time.time() - t0
print(f"  Searched in {elapsed:.1f}s, found {len(best)} candidates")

best.sort(key=lambda x: -x[0])
for C, Cd, coeffs in best[:20]:
    poly_str = f"x⁵"
    if coeffs[1]: poly_str += f" + {coeffs[1]}x⁴" if coeffs[1] > 0 else f" - {-coeffs[1]}x⁴"
    if coeffs[2]: poly_str += f" + {coeffs[2]}x³" if coeffs[2] > 0 else f" - {-coeffs[2]}x³"
    if coeffs[3]: poly_str += f" + {coeffs[3]}x²" if coeffs[3] > 0 else f" - {-coeffs[3]}x²"
    if coeffs[4]: poly_str += f" + {coeffs[4]}x" if coeffs[4] > 0 else f" - {-coeffs[4]}x"
    if coeffs[5]: poly_str += f" + {coeffs[5]}" if coeffs[5] > 0 else f" - {-coeffs[5]}"
    # Check irreducibility
    sym_poly = sum(c * x**i for i, c in enumerate(reversed(coeffs)))
    fac = factor(sym_poly, domain=QQ)
    irr = "IRR" if str(fac) == str(sym_poly) else "red"
    print(f"  C_f={C:>7.3f}  C/d={Cd:>6.3f}  {irr}  {poly_str}")


# ── Strategy 3: Smarter construction ─────────────────────────────
# Use discriminant theory: for a quintic to block prime q,
# we need f(x) to have no roots mod q.
# Maximise blocked primes by choosing f so its discriminant
# has special properties.

print(f"\n{'='*70}")
print("SEARCH 2: QUINTICS WITH LEADING COEFF > 1")
print(f"{'='*70}")

best2 = []
t0 = time.time()

for lead in [2, 3, 4, 5, 6, 8]:
    for a in range(-5, 6):
        for b in range(-15, 16):
            for c in range(-30, 31, 2):
                for d in range(-40, 41, 3):
                    for e in [1, -1, 3, -3, 7, -7, 11, -11, 41, -41]:
                        coeffs = [lead, a, b, c, d, e]
                        blocked = sum(1 for p in screen_primes[:15]
                                      if count_roots_mod_p(coeffs, p) == 0)
                        if blocked < 11:
                            continue
                        C = bateman_horn_fast(coeffs, screen_primes)
                        if C / 5 > euler_Cd * 0.8:
                            C_full = bateman_horn_fast(coeffs)
                            if C_full / 5 > euler_Cd * 0.7:
                                best2.append((C_full, C_full/5, coeffs))

elapsed = time.time() - t0
print(f"  Searched in {elapsed:.1f}s, found {len(best2)} candidates")

best2.sort(key=lambda x: -x[0])
for C, Cd, coeffs in best2[:20]:
    sym_poly = sum(c * x**i for i, c in enumerate(reversed(coeffs)))
    fac = factor(sym_poly, domain=QQ)
    irr = "IRR" if str(fac) == str(sym_poly) else "red"
    print(f"  C_f={C:>7.3f}  C/d={Cd:>6.3f}  {irr}  {sym_poly}")

# ── Evaluate the best candidates ─────────────────────────────────

print(f"\n{'='*70}")
print("TOP CANDIDATES: ACTUAL PRIME DENSITY AT PRIME ARGUMENTS")
print(f"{'='*70}")

all_best = best + best2
all_best.sort(key=lambda x: -x[0])

# Deduplicate
seen = set()
unique_best = []
for C, Cd, coeffs in all_best:
    key = tuple(coeffs)
    if key not in seen:
        seen.add(key)
        unique_best.append((C, Cd, coeffs))

print(f"\n  {'C_f':>7} {'C/d':>6} {'primes':>7} {'density':>8} {'irr?':>5}  Polynomial")
print(f"  {'-'*70}")

# Always show Euler for comparison
euler_count, euler_dens = prime_density([1, 1, 41])
print(f"  {euler_C:>7.3f} {euler_Cd:>6.3f} {euler_count:>7} {euler_dens:>8.3f} {'IRR':>5}  p²+p+41 (EULER)")

for C, Cd, coeffs in unique_best[:15]:
    sym_poly = sum(c * x**i for i, c in enumerate(reversed(coeffs)))
    fac = factor(sym_poly, domain=QQ)
    is_irr = str(fac) == str(sym_poly)
    if not is_irr:
        continue  # skip reducible
    count, dens = prime_density(coeffs)
    beat = " ★ BEATS EULER" if dens > euler_dens else ""
    print(f"  {C:>7.3f} {Cd:>6.3f} {count:>7} {dens:>8.3f} {'IRR':>5}  {sym_poly}{beat}")
