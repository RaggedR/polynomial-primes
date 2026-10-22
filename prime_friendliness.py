"""
How "prime-friendly" is each irreducible Chebyshev polynomial T_{2^k}?

Compute:
1. The Bateman-Horn constant C_f (Euler product over small primes)
2. Actual prime counts for T_{2^k}(n), n = 1..N
3. Comparison of predicted vs actual
4. Which small primes are blocked by quadratic reciprocity generalisation?
"""

import math
from sympy import chebyshevt, Symbol, Poly, QQ, isprime, primerange

x = Symbol('x')

# Precompute polynomial expressions
polys = {}
for k in range(5):  # T_1, T_2, T_4, T_8, T_16
    n = 2**k
    T = chebyshevt(n, x)
    polys[n] = T
    print(f"T_{n}(x) = {T}")

print()

# ── Bateman-Horn constant ────────────────────────────────────────

def count_roots_mod_p(f_poly, p):
    """Count solutions of f(x) ≡ 0 (mod p)."""
    count = 0
    coeffs = Poly(f_poly, x, domain=QQ).all_coeffs()
    for x_val in range(p):
        val = 0
        for c in coeffs:
            val = (val * x_val + int(c)) % p
        if val == 0:
            count += 1
    return count


print("=" * 72)
print("BATEMAN-HORN CONSTANTS")
print("=" * 72)

sieve_primes = list(primerange(2, 200))

for k in range(5):
    n = 2**k
    T = polys[n]
    deg = 2**k

    print(f"\n  T_{n}(x), degree {deg}")
    print(f"  {'p':>5} {'ω(p)':>5} {'factor':>10}")
    print(f"  {'-'*25}")

    C = 1.0
    blocked = []
    for p in sieve_primes:
        w = count_roots_mod_p(T, p)
        if p == deg * (2 ** (k+1)):  # leading coeff might cause issues
            pass
        factor = (1 - w / p) / (1 - 1 / p)
        C *= factor
        if p < 30:
            blocked_str = " ← BLOCKED" if w == 0 else ""
            print(f"  {p:>5} {w:>5} {factor:>10.4f}{blocked_str}")
        if w == 0 and p < 50:
            blocked.append(p)

    print(f"\n  C_f ≈ {C:.4f}")
    print(f"  Primes < 50 that can NEVER divide T_{n}(m): {blocked}")


# ── Actual prime counts ──────────────────────────────────────────

print("\n" + "=" * 72)
print("ACTUAL PRIME COUNTS: T_{2^k}(n) for n = 1..N")
print("=" * 72)

def eval_T(poly_expr, val):
    """Evaluate polynomial at integer value."""
    return int(poly_expr.subs(x, val))

for k in range(5):
    n = 2**k
    T = polys[n]
    deg = n

    print(f"\n  T_{n}(x), degree {deg}:")

    # Compute C_f for prediction
    C = 1.0
    for p in sieve_primes:
        w = count_roots_mod_p(T, p)
        factor = (1 - w / p) / (1 - 1 / p)
        C *= factor

    print(f"  C_f ≈ {C:.4f}")
    print(f"  {'N':>8} {'predicted':>10} {'actual':>8} {'ratio':>8} {'density':>10}")

    for N in [50, 100, 500, 1000, 5000]:
        predicted = C * N / (deg * math.log(N))
        count = 0
        for m in range(2, N + 1):  # skip m=0,1 for cleanliness
            val = eval_T(T, m)
            if abs(val) > 1 and isprime(abs(val)):
                count += 1
        density = count / N
        ratio = count / predicted if predicted > 0 else 0
        print(f"  {N:>8} {predicted:>10.1f} {count:>8} {ratio:>8.3f} {density:>10.3f}")


# ── First prime values of each ───────────────────────────────────

print("\n" + "=" * 72)
print("FIRST PRIME VALUES")
print("=" * 72)

for k in range(5):
    n = 2**k
    T = polys[n]
    print(f"\n  T_{n}(x):")

    found = 0
    for m in range(2, 200):
        val = eval_T(T, m)
        if abs(val) > 1 and isprime(abs(val)):
            print(f"    T_{n}({m}) = {val}")
            found += 1
            if found >= 8:
                break


# ── The "shield" visualisation ───────────────────────────────────

print("\n" + "=" * 72)
print("SMALL PRIME SHIELD: which primes q can divide T_{2^k}(n)?")
print("=" * 72)
print(f"\n  {'q':>4}", end="")
for k in range(5):
    n = 2**k
    print(f"{'T_'+str(n):>8}", end="")
print()
print(f"  {'':>4}", end="")
for k in range(5):
    print(f"{'ω(q)':>8}", end="")
print()
print(f"  {'-'*44}")

for q in primerange(2, 60):
    print(f"  {q:>4}", end="")
    for k in range(5):
        n = 2**k
        T = polys[n]
        w = count_roots_mod_p(T, q)
        if w == 0:
            print(f"{'·':>8}", end="")  # blocked
        else:
            print(f"{w:>8}", end="")
    print()

print(f"\n  · = blocked (ω = 0, this prime can never divide T_n(m))")
print(f"  Numbers = how many residue classes mod q allow divisibility")


# ── Summary ──────────────────────────────────────────────────────

print("\n" + "=" * 72)
print("SUMMARY: PRIME-FRIENDLINESS RANKING")
print("=" * 72)

results = []
for k in range(5):
    n = 2**k
    T = polys[n]
    deg = n
    C = 1.0
    for p in sieve_primes:
        w = count_roots_mod_p(T, p)
        C *= (1 - w / p) / (1 - 1 / p)

    blocked_count = sum(1 for p in primerange(2, 50) if count_roots_mod_p(T, p) == 0)
    results.append((n, deg, C, blocked_count))

print(f"\n  {'T_n':>5} {'degree':>7} {'C_f':>8} {'blocked<50':>12} {'C_f/degree':>11}")
print(f"  {'-'*50}")
for n, deg, C, bc in results:
    # C_f / degree is the effective prime generation rate
    print(f"  {'T_'+str(n):>5} {deg:>7} {C:>8.3f} {bc:>12} {C/deg:>11.4f}")

print(f"""
  C_f/degree is the effective rate: predicted primes ~ (C_f/degree) · N/ln(N)
  Higher = more prime-friendly.
""")
