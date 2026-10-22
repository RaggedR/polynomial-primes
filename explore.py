"""
Exploring classical polynomial families evaluated at primes.

Questions:
1. When is T_n(p) prime? (Chebyshev of first kind)
2. What divisibility patterns emerge in T_n(p), H_n(p)?
3. Does gcd(T_m(p), T_n(p)) relate to T_{gcd(m,n)}(p)?
4. Do Hermite values H_n(p) have interesting prime factorisations?
"""

from sympy import (
    isprime, factorint, gcd, primerange, chebyshevt, hermite,
    legendre, Symbol, Rational
)
from math import gcd as mgcd
from collections import defaultdict

x = Symbol('x')

# Precompute polynomial functions
def T(n, val):
    """Chebyshev of the first kind T_n(val)."""
    return int(chebyshevt(n, x).subs(x, val))

def H(n, val):
    """Physicist's Hermite H_n(val)."""
    return int(hermite(n, x).subs(x, val))

primes = list(primerange(2, 50))

# ── 1. Chebyshev T_n(p): when is the result prime? ──────────────

print("=" * 70)
print("CHEBYSHEV T_n(p): PRIMALITY")
print("=" * 70)

for n in range(2, 8):
    poly_expr = chebyshevt(n, x)
    print(f"\nT_{n}(x) = {poly_expr}")
    results = []
    for p in primes:
        val = T(n, p)
        if isprime(abs(val)):
            results.append((p, val))
    if results:
        print(f"  Prime values at primes: {results[:10]}")
    else:
        print(f"  No prime values at primes in range")

# ── 2. Structural divisibility of T_n(p) ────────────────────────

print("\n" + "=" * 70)
print("CHEBYSHEV T_n(p): DIVISIBILITY BY p")
print("=" * 70)

for n in range(2, 10):
    divides = all(T(n, p) % p == 0 for p in primes)
    print(f"  p | T_{n}(p) for all primes p tested: {divides}")

# ── 3. T_2(p) = 2p² - 1: a prime-generating machine? ───────────

print("\n" + "=" * 70)
print("T_2(p) = 2p² - 1: PRIME DENSITY")
print("=" * 70)

bigger_primes = list(primerange(2, 200))
t2_primes = []
t2_composites = []
for p in bigger_primes:
    val = 2 * p * p - 1
    if isprime(val):
        t2_primes.append((p, val))
    else:
        t2_composites.append((p, val, factorint(val)))

print(f"  Primes p < 200 where T_2(p) is prime: {len(t2_primes)}/{len(bigger_primes)}")
print(f"  Density: {len(t2_primes)/len(bigger_primes):.3f}")
print(f"  First few: {t2_primes[:15]}")
print(f"\n  First few composite T_2(p):")
for p, val, factors in t2_composites[:8]:
    print(f"    T_2({p}) = {val} = {factors}")

# ── 4. GCD structure: gcd(T_m(p), T_n(p)) vs T_{gcd(m,n)}(p) ──

print("\n" + "=" * 70)
print("GCD STRUCTURE: gcd(T_m(p), T_n(p)) vs T_{gcd(m,n)}(p)")
print("=" * 70)

for p in [2, 3, 5, 7]:
    print(f"\n  p = {p}:")
    for m in range(2, 7):
        for n in range(m + 1, 7):
            g = gcd(T(m, p), T(n, p))
            t_gcd = T(mgcd(m, n), p)
            match = "✓" if g == abs(t_gcd) else f"✗ (gcd={g}, T_{mgcd(m,n)}={t_gcd})"
            if g != abs(t_gcd):
                # Check if T_{gcd} divides the gcd
                divides = g % abs(t_gcd) == 0 if t_gcd != 0 else "N/A"
                print(f"    gcd(T_{m}, T_{n}) = {g}, T_{mgcd(m,n)}({p}) = {t_gcd}, "
                      f"T_gcd | gcd: {divides}")

# ── 5. Hermite H_n(p): factorisation patterns ──────────────────

print("\n" + "=" * 70)
print("HERMITE H_n(p): FACTORISATION PATTERNS")
print("=" * 70)

for n in range(2, 7):
    poly_expr = hermite(n, x)
    print(f"\nH_{n}(x) = {poly_expr}")
    for p in primes[:6]:
        val = H(n, p)
        if val != 0:
            factors = factorint(abs(val))
            sign = "-" if val < 0 else ""
            print(f"  H_{n}({p}) = {sign}{factors}")

# ── 6. v_p(H_n(p)): p-adic valuation of Hermite at p ───────────

print("\n" + "=" * 70)
print("p-ADIC VALUATION: v_p(H_n(p))")
print("=" * 70)

def v_p(n, p):
    """p-adic valuation of n."""
    if n == 0:
        return float('inf')
    n = abs(n)
    count = 0
    while n % p == 0:
        n //= p
        count += 1
    return count

for p in [2, 3, 5, 7]:
    vals = []
    for n in range(1, 16):
        val = H(n, p)
        vals.append(v_p(val, p))
    print(f"  p={p}: v_p(H_n(p)) for n=1..15: {vals}")

# ── 7. Chebyshev T_n(p) mod small primes ────────────────────────

print("\n" + "=" * 70)
print("T_n(p) mod q: RESIDUE PATTERNS")
print("=" * 70)
print("(rows = n, columns = p, entries = T_n(p) mod 7)")

print(f"\n  {'n\\p':>4}", end="")
for p in primes[:10]:
    print(f"{p:>5}", end="")
print()

for n in range(1, 10):
    print(f"  {n:>4}", end="")
    for p in primes[:10]:
        print(f"{T(n, p) % 7:>5}", end="")
    print()

# ── 8. Connection to Lucas sequences ────────────────────────────

print("\n" + "=" * 70)
print("LUCAS CONNECTION: 2·T_n(p/2) vs Lucas V_n(p, 1)")
print("=" * 70)
print("Lucas V_n(P, Q) with P=p, Q=1 satisfies V_n = P·V_{n-1} - Q·V_{n-2}")

def lucas_V(n, P, Q=1):
    """Lucas sequence V_n(P, Q)."""
    if n == 0:
        return 2
    if n == 1:
        return P
    v_prev, v_curr = 2, P
    for _ in range(2, n + 1):
        v_prev, v_curr = v_curr, P * v_curr - Q * v_prev
    return v_curr

for p in [2, 3, 5, 7, 11]:
    print(f"\n  p = {p}:")
    for n in range(1, 8):
        t_val = T(n, Rational(p, 2))
        lucas_val = lucas_V(n, p)
        scaled_t = 2 * t_val
        match = "✓" if scaled_t == lucas_val else "✗"
        print(f"    2·T_{n}({p}/2) = {scaled_t}, V_{n}({p},1) = {lucas_val}  {match}")

# ── 9. T_n(p) and quadratic residues ────────────────────────────

print("\n" + "=" * 70)
print("QUADRATIC RESIDUES: is T_n(p) a QR mod q?")
print("=" * 70)
from sympy.ntheory import jacobi_symbol

for q in [3, 5, 7, 11, 13]:
    print(f"\n  Jacobi symbol (T_n(p) / {q}) for p=2, n=1..12:")
    symbols = []
    for n in range(1, 13):
        val = T(n, 2)
        if val % q != 0:
            symbols.append(jacobi_symbol(val, q))
        else:
            symbols.append(0)
    print(f"    {symbols}")
    # Check periodicity
    for period in range(1, 7):
        if all(symbols[i] == symbols[i + period] for i in range(len(symbols) - period)):
            print(f"    Period {period} detected!")
            break

print("\n" + "=" * 70)
print("SUMMARY OF FINDINGS")
print("=" * 70)
print("""
Key observations to investigate:
1. T_2(p) = 2p²-1 generates primes at high density
2. For odd n, p | T_n(p) always (structural)
3. 2·T_n(p/2) = V_n(p,1) — Chebyshev IS the Lucas sequence
4. The GCD structure and periodicity mod q
5. p-adic valuations of H_n(p) show regular patterns
""")
