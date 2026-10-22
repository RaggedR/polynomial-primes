"""
"Chebyshev primes": values of T_2(p) = 2p² - 1 for prime p.
Question: does every positive integer have a unique decomposition
as a sum of Chebyshev primes (with repetition allowed)?

Also explore: what if we only use the ones that are actually prime?
"""

from sympy import primerange, isprime
from collections import defaultdict

# ── Build the set of "Chebyshev primes" ──────────────────────────

primes_list = list(primerange(2, 100))

# All T_2(p) values
cheb_all = sorted(set(2 * p * p - 1 for p in primes_list))
# Only the ones that are actually prime
cheb_prime = sorted(c for c in cheb_all if isprime(c))

print("Chebyshev numbers (all T_2(p)):")
print(f"  {cheb_all[:15]}...")
print(f"\nChebyshev primes (T_2(p) where result is prime):")
print(f"  {cheb_prime[:15]}...")

# ── Count representations: how many ways to write n as a sum? ────

def count_representations(target, basis, max_terms=None):
    """
    Count the number of ways to write target as an ordered sum
    of elements from basis (with repetition). For unique decomposition,
    we want UNordered (partition-style).

    Returns dict: n -> number of unordered representations.
    Uses dynamic programming.
    """
    # Number of ways to write each n as unordered sum of elements from basis
    ways = [0] * (target + 1)
    ways[0] = 1  # empty sum

    # Unordered: process each basis element in turn (like coin change)
    for b in basis:
        for n in range(b, target + 1):
            ways[n] += ways[n - b]

    return ways


N = 2000

print(f"\n{'=' * 70}")
print(f"ADDITIVE DECOMPOSITION INTO CHEBYSHEV PRIMES (≤ {N})")
print(f"{'=' * 70}")

# Use only Chebyshev primes that are ≤ N
basis_prime = [c for c in cheb_prime if c <= N]
basis_all = [c for c in cheb_all if c <= N]

print(f"\nBasis (primes only): {basis_prime}")
print(f"Basis (all T_2(p)):  {basis_all}")

for label, basis in [("CHEBYSHEV PRIMES (prime T_2(p) only)", basis_prime),
                     ("ALL CHEBYSHEV NUMBERS (all T_2(p))", basis_all)]:
    ways = count_representations(N, basis)

    representable = sum(1 for n in range(1, N+1) if ways[n] > 0)
    unique = sum(1 for n in range(1, N+1) if ways[n] == 1)
    multi = sum(1 for n in range(1, N+1) if ways[n] > 1)
    unrep = sum(1 for n in range(1, N+1) if ways[n] == 0)

    print(f"\n  {label}:")
    print(f"  Representable:     {representable}/{N} ({representable/N:.1%})")
    print(f"  Unique decomp:     {unique}/{N}")
    print(f"  Multiple decomps:  {multi}/{N}")
    print(f"  Not representable: {unrep}/{N}")

    # Show first few unrepresentable
    unrep_list = [n for n in range(1, N+1) if ways[n] == 0]
    if unrep_list:
        print(f"  First unrepresentable: {unrep_list[:20]}")

    # Show first few with unique decomposition
    unique_list = [n for n in range(1, N+1) if ways[n] == 1]
    if unique_list:
        print(f"  First with unique decomp: {unique_list[:20]}")

    # Show numbers with most representations
    max_ways = max(ways[n] for n in range(1, N+1))
    most_reps = [(n, ways[n]) for n in range(1, N+1) if ways[n] == max_ways]
    print(f"  Max representations: {max_ways} (at n={most_reps[0][0]})")

    # Distribution of representation counts
    dist = defaultdict(int)
    for n in range(1, N+1):
        w = ways[n]
        if w == 0:
            dist[0] += 1
        elif w == 1:
            dist[1] += 1
        elif w <= 5:
            dist['2-5'] += 1
        elif w <= 20:
            dist['6-20'] += 1
        else:
            dist['21+'] += 1

    print(f"  Distribution: 0 reps: {dist[0]}, 1 rep: {dist[1]}, "
          f"2-5: {dist.get('2-5',0)}, 6-20: {dist.get('6-20',0)}, "
          f"21+: {dist.get('21+',0)}")


# ── Find actual decompositions for small numbers ─────────────────

def find_decomposition(n, basis):
    """Find one unordered decomposition of n into elements of basis."""
    if n == 0:
        return []
    if n < 0:
        return None
    for b in reversed(basis):  # greedy: try largest first
        if b <= n:
            rest = find_decomposition(n - b, [x for x in basis if x <= b])
            if rest is not None:
                return [b] + rest
    return None

def find_all_decompositions(n, basis, max_elem=None):
    """Find all unordered decompositions."""
    if max_elem is None:
        max_elem = max(basis) if basis else 0
    if n == 0:
        return [[]]
    if n < 0:
        return []
    results = []
    for b in basis:
        if b <= n and b <= max_elem:
            for rest in find_all_decompositions(n - b, basis, b):
                results.append([b] + rest)
    return results


print(f"\n{'=' * 70}")
print(f"DECOMPOSITIONS OF SMALL NUMBERS INTO CHEBYSHEV PRIMES")
print(f"{'=' * 70}")

basis = basis_prime
for n in range(1, 80):
    w = count_representations(80, basis)[n]
    if w == 0:
        continue
    decomps = find_all_decompositions(n, basis)
    if len(decomps) == 1:
        d = decomps[0]
        print(f"  {n:>3} = {' + '.join(str(x) for x in d)}  (unique ✓)")
    else:
        print(f"  {n:>3} = {len(decomps)} ways:")
        for d in decomps[:4]:
            print(f"        {' + '.join(str(x) for x in d)}")
        if len(decomps) > 4:
            print(f"        ... and {len(decomps)-4} more")


# ── What about multiplicative decomposition? ─────────────────────

print(f"\n{'=' * 70}")
print(f"MULTIPLICATIVE STRUCTURE")
print(f"{'=' * 70}")
print("Can every integer be written as a PRODUCT of Chebyshev primes?")
print("(This would be more analogous to unique prime factorisation)")

def mult_decompose(n, basis, max_elem=None):
    """Find all unordered factorisations of n using elements of basis."""
    if max_elem is None:
        max_elem = max(basis) if basis else 0
    if n == 1:
        return [[]]
    results = []
    for b in basis:
        if b <= max_elem and b <= n and n % b == 0:
            for rest in mult_decompose(n // b, basis, b):
                results.append([b] + rest)
    return results

print(f"\nMultiplicative decompositions (using Chebyshev primes {basis_prime[:8]}):")
has_decomp = 0
unique_decomp = 0
for n in range(2, 500):
    decomps = mult_decompose(n, basis_prime)
    if decomps:
        has_decomp += 1
        if len(decomps) == 1:
            unique_decomp += 1
        if n <= 120 or len(decomps) > 1:
            if decomps and n <= 120:
                for d in decomps:
                    print(f"  {n:>5} = {' × '.join(str(x) for x in d)}"
                          f"{'  (unique ✓)' if len(decomps) == 1 else ''}")

print(f"\n  Of integers 2..499: {has_decomp} have a Chebyshev factorisation, "
      f"{unique_decomp} unique")


# ── Greedy representation (like binary/Zeckendorf) ───────────────

print(f"\n{'=' * 70}")
print(f"GREEDY REPRESENTATION (largest-first, no repeats)")
print(f"{'=' * 70}")
print("Like Zeckendorf (Fibonacci) — use each Chebyshev prime at most once.")
print("Is this always possible? Is it unique?\n")

def greedy_no_repeat(n, basis):
    """Greedy: use each element at most once, largest first."""
    used = []
    remaining = n
    for b in reversed(basis):
        if b <= remaining:
            used.append(b)
            remaining -= b
    if remaining == 0:
        return used
    return None  # greedy failed

def count_subset_sums(target, basis):
    """Count subsets of basis that sum to target (each used at most once)."""
    ways = [0] * (target + 1)
    ways[0] = 1
    for b in basis:
        for n in range(target, b - 1, -1):
            ways[n] += ways[n - b]
    return ways

basis_big = [c for c in cheb_prime if c <= 5000]
ways_no_rep = count_subset_sums(2000, basis_big)

representable = sum(1 for n in range(1, 2001) if ways_no_rep[n] > 0)
unique = sum(1 for n in range(1, 2001) if ways_no_rep[n] == 1)

print(f"Using Chebyshev primes ≤ 5000, no repeats, for n ≤ 2000:")
print(f"  Representable: {representable}/2000")
print(f"  Unique:        {unique}/2000")

# Show some greedy decompositions
for n in [42, 100, 200, 500, 1000]:
    g = greedy_no_repeat(n, basis_big)
    w = ways_no_rep[n] if n <= 2000 else "?"
    if g:
        print(f"\n  {n} = {' + '.join(str(x) for x in g)}  ({w} total ways)")
    else:
        print(f"\n  {n}: greedy fails  ({w} total ways)")
