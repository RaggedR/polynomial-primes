"""
Can we use 2p² - 1 as a practical large prime generator?

Strategy: pick successive large primes p, compute 2p²-1, test primality.
Compare: how many candidates do we test vs random odd numbers of the same size?
"""

import time
from sympy import isprime, nextprime, primerange
import random
import math

# ── Method 1: Chebyshev T_2(p) = 2p² - 1 ────────────────────────

def generate_prime_chebyshev(start_p):
    """Generate a prime of the form 2p²-1, starting search at prime p."""
    p = nextprime(start_p - 1)
    attempts = 0
    while True:
        candidate = 2 * p * p - 1
        attempts += 1
        if isprime(candidate):
            return candidate, p, attempts
        p = nextprime(p)


# ── Method 2: Random odd numbers of similar size ─────────────────

def generate_prime_random(target_size):
    """Generate a random prime near target_size."""
    attempts = 0
    while True:
        candidate = random.randrange(target_size // 2, target_size * 2) | 1
        attempts += 1
        if isprime(candidate):
            return candidate, attempts


# ── Method 3: Meixner-Pollaczek 8p² - 31 ────────────────────────

def generate_prime_mp(start_p):
    """Generate a prime of the form 8p²-31."""
    p = nextprime(start_p - 1)
    attempts = 0
    while True:
        candidate = 8 * p * p - 31
        attempts += 1
        if candidate > 1 and isprime(candidate):
            return candidate, p, attempts
        p = nextprime(p)


# ── Benchmark at various sizes ───────────────────────────────────

print("=" * 75)
print("LARGE PRIME GENERATION: CHEBYSHEV vs MEIXNER-POLLACZEK vs RANDOM")
print("=" * 75)

print(f"\n{'start p':>12} {'digits':>7} {'method':>8} {'attempts':>9} {'time (s)':>9} {'prime found'}")
print("-" * 75)

for start in [10**3, 10**4, 10**5, 10**6, 10**7]:
    target_size = 2 * start * start  # approximate size of 2p²-1

    # Chebyshev
    t0 = time.time()
    prime_found, p_used, attempts = generate_prime_chebyshev(start)
    t1 = time.time()
    digits = len(str(prime_found))
    print(f"{start:>12} {digits:>7} {'T_2':>8} {attempts:>9} {t1-t0:>9.3f} {str(prime_found)[:30]}...")

    # Meixner-Pollaczek
    t0 = time.time()
    prime_found_mp, p_used_mp, attempts_mp = generate_prime_mp(start)
    t1 = time.time()
    digits_mp = len(str(prime_found_mp))
    print(f"{'':>12} {digits_mp:>7} {'MP_7':>8} {attempts_mp:>9} {t1-t0:>9.3f} {str(prime_found_mp)[:30]}...")

    # Random
    t0 = time.time()
    prime_found_r, attempts_r = generate_prime_random(target_size)
    t1 = time.time()
    digits_r = len(str(prime_found_r))
    print(f"{'':>12} {digits_r:>7} {'random':>8} {attempts_r:>9} {t1-t0:>9.3f} {str(prime_found_r)[:30]}...")

    # Expected attempts
    expected_random = math.log(target_size) / 2  # ~ln(N)/2 for odd numbers
    print(f"{'':>12} {'':>7} {'expect':>8} {'~'+str(int(expected_random)):>9}")
    print()


# ── Theoretical comparison ───────────────────────────────────────

print("=" * 75)
print("WHY THIS WORKS (AND WHY IT DOESN'T)")
print("=" * 75)

print("""
For a number near N = 2p²:
  Random odd number:  Pr[prime] ≈ 2/ln(N) ≈ 1/ln(p)
  2p²-1 at prime p:   Pr[prime] ≈ C_f/(2·ln(p))  with C_f ≈ 3.68

So the Chebyshev method needs ~2·ln(p)/3.68 attempts,
vs random needing ~ln(N)/2 = ln(p) + ln(2)/2 ≈ ln(p) attempts.

Ratio: Chebyshev needs ~55% as many attempts as random.
       MP_7 needs ~52% as many attempts (C_f ≈ 3.88).

BUT: each attempt costs the SAME. Primality testing of 2p²-1
uses Baillie-PSW just like any other number. There's no shortcut
like Lucas-Lehmer for Mersenne numbers.
""")

# ── Compare with Mersenne ────────────────────────────────────────

print("=" * 75)
print("THE MERSENNE ADVANTAGE")
print("=" * 75)

print("""
Mersenne primes M_p = 2^p - 1 have a STRUCTURAL advantage:
  - Lucas-Lehmer test: O(p² log p) — tests M_p in time proportional to p²
  - General primality test: O(N^(1/3+ε)) at best for a number with N digits

For a 1000-digit number:
  - Mersenne LL test:  ~3300² ≈ 10^7 operations
  - Baillie-PSW:       ~1000³ ≈ 10^9 operations (modular exponentiations)

So Mersenne is ~100x faster PER TEST, but primes are exponentially rarer
(density ~1/p vs ~1/ln(p) for quadratics).

The Chebyshev method's advantage is DENSITY, not speed:
  - ~1 in 5 candidates is prime (at moderate sizes)
  - vs ~1 in ln(N)/2 for random numbers
  - but the per-test cost is identical
""")

# ── Can we do better? Use 2p²-1 where p is ALSO of this form? ──

print("=" * 75)
print("ITERATED CHEBYSHEV: TOWERS OF PRIMES")
print("=" * 75)
print("Start with a small prime, iterate p → 2p²-1. How far can we go?\n")

p = 2
chain = [p]
for i in range(20):
    candidate = 2 * p * p - 1
    if isprime(candidate):
        p = candidate
        chain.append(p)
        print(f"  Step {i+1}: 2·{chain[-2]}² - 1 = {p}  ✓ PRIME  ({len(str(p))} digits)")
    else:
        print(f"  Step {i+1}: 2·{p}² - 1 = {candidate}  ✗ composite")
        # Try next prime after p
        break

print(f"\n  Chain died at length {len(chain)}: {chain}")

# Try multiple starting primes
print("\nChains from different starting primes:")
for start in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
    p = start
    chain = [p]
    for i in range(15):
        candidate = 2 * p * p - 1
        if candidate > 10**30:
            break
        if isprime(candidate):
            p = candidate
            chain.append(p)
        else:
            break
    print(f"  Start {start:>3}: chain length {len(chain):>2}  →  {chain}")
