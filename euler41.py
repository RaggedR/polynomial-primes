"""
Does Euler's n² + n + 41 beat our winner 6n² - 6n + 1?
"""

from sympy import isprime, primerange
import math

primes_list = list(primerange(2, 10000))

polys = [
    ("6p²-6p+1  (Bernoulli B₂)",   lambda p: 6*p*p - 6*p + 1),
    ("p²+p+41   (Euler)",           lambda p: p*p + p + 41),
    ("2p²-1     (Chebyshev T₂)",    lambda p: 2*p*p - 1),
    ("p²+1",                        lambda p: p*p + 1),
    ("p²+p+17   (Euler's other)",   lambda p: p*p + p + 17),
]

# ── Bateman-Horn constants ───────────────────────────────────────

sieve_primes = list(primerange(2, 500))

def count_roots(f_coeffs, p):
    """f_coeffs = [a, b, c] for ax²+bx+c. Count roots mod p."""
    a, b, c = f_coeffs
    count = 0
    for x in range(p):
        if (a*x*x + b*x + c) % p == 0:
            count += 1
    return count

def bateman_horn(coeffs):
    C = 1.0
    for p in sieve_primes:
        w = count_roots(coeffs, p)
        C *= (1 - w/p) / (1 - 1/p)
    return C

poly_coeffs = [
    ("6p²-6p+1  (Bernoulli B₂)",  (6, -6, 1)),
    ("p²+p+41   (Euler)",          (1, 1, 41)),
    ("2p²-1     (Chebyshev T₂)",   (2, 0, -1)),
    ("p²+1",                       (1, 0, 1)),
    ("p²+p+17   (Euler's other)",  (1, 1, 17)),
]

print("=" * 75)
print("BATEMAN-HORN CONSTANTS AND DISCRIMINANTS")
print("=" * 75)
print(f"\n  {'Polynomial':>30} {'Disc':>6} {'C_f':>7} {'C_f/2':>7}")
print(f"  {'-'*55}")

for name, (a, b, c) in poly_coeffs:
    disc = b*b - 4*a*c
    C = bateman_horn((a, b, c))
    print(f"  {name:>30} {disc:>6} {C:>7.3f} {C/2:>7.4f}")

# ── Blocked primes ───────────────────────────────────────────────

print(f"\n{'=' * 75}")
print("SHIELD: which primes < 50 are blocked?")
print("=" * 75)

for name, (a, b, c) in poly_coeffs:
    blocked = []
    for p in primerange(2, 50):
        if count_roots((a, b, c), p) == 0:
            blocked.append(p)
    print(f"  {name:>30}: {blocked}")

# ── Prime counts at prime arguments ──────────────────────────────

print(f"\n{'=' * 75}")
print("PRIME COUNTS AT PRIME ARGUMENTS p < 10000")
print("=" * 75)

print(f"\n  {'Polynomial':>30} {'primes':>7} {'/ 1229':>7} {'density':>8}")
print(f"  {'-'*55}")

for name, f in polys:
    count = sum(1 for p in primes_list if abs(f(p)) > 1 and isprime(abs(f(p))))
    print(f"  {name:>30} {count:>7} {'/ 1229':>7} {count/len(primes_list):>8.3f}")

# ── Head to head table ───────────────────────────────────────────

print(f"\n{'=' * 75}")
print("HEAD TO HEAD: first 20 primes")
print("=" * 75)

print(f"\n  {'p':>4} {'6p²-6p+1':>12} {'?':>2} {'p²+p+41':>12} {'?':>2} {'2p²-1':>10} {'?':>2}")
print(f"  {'-'*50}")

for p in list(primerange(2, 75)):
    v1 = 6*p*p - 6*p + 1
    v2 = p*p + p + 41
    v3 = 2*p*p - 1
    p1 = "✓" if isprime(v1) else ""
    p2 = "✓" if isprime(v2) else ""
    p3 = "✓" if isprime(v3) else ""
    print(f"  {p:>4} {v1:>12} {p1:>2} {v2:>12} {p2:>2} {v3:>10} {p3:>2}")
