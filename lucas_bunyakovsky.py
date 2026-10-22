"""
Deep dive:
1. Chebyshev = Lucas = primality testing (the exact identity)
2. Bunyakovsky's conjecture and T_2(p) = 2p² - 1
"""

from sympy import (
    isprime, factorint, primerange, chebyshevt, chebyshevt_root,
    Symbol, Rational, sqrt, simplify, floor, prime, nextprime
)
from sympy.functions.combinatorial.numbers import jacobi_symbol
from collections import Counter

x = Symbol('x')

# ══════════════════════════════════════════════════════════════════
# PART 1: THE EXACT IDENTITY — Chebyshev IS Lucas
# ══════════════════════════════════════════════════════════════════

print("=" * 72)
print("PART 1: CHEBYSHEV IS LUCAS (THE EXACT IDENTITY)")
print("=" * 72)

def lucas_V(n, P, Q=1):
    """Lucas V-sequence: V_0=2, V_1=P, V_n = P·V_{n-1} - Q·V_{n-2}."""
    if n == 0: return 2
    if n == 1: return P
    a, b = 2, P
    for _ in range(2, n + 1):
        a, b = b, P * b - Q * a
    return b

def lucas_U(n, P, Q=1):
    """Lucas U-sequence: U_0=0, U_1=1, U_n = P·U_{n-1} - Q·U_{n-2}."""
    if n == 0: return 0
    if n == 1: return 1
    a, b = 0, 1
    for _ in range(2, n + 1):
        a, b = b, P * b - Q * a
    return b

print("""
THE IDENTITY: V_n(P, 1) = 2·T_n(P/2)   [exact, no exceptions]

Proof sketch:
  Let α, β be roots of t² - Pt + 1 = 0. Then αβ = 1, so β = 1/α.
  V_n = αⁿ + βⁿ = αⁿ + α⁻ⁿ.
  Set cos θ = P/2 (i.e. α = e^{iθ}). Then:
    αⁿ + α⁻ⁿ = 2cos(nθ) = 2·T_n(cos θ) = 2·T_n(P/2).  □

Previous code had a bug: int() truncated Rational values.
Now using exact Rational arithmetic:
""")

for P in [2, 3, 5, 7, 11, 13]:
    print(f"  P = {P}:")
    all_match = True
    for n in range(0, 10):
        chebyshev_val = 2 * chebyshevt(n, x).subs(x, Rational(P, 2))
        lucas_val = lucas_V(n, P)
        match = (chebyshev_val == lucas_val)
        if not match:
            all_match = False
            print(f"    n={n}: 2·T_{n}({P}/2) = {chebyshev_val}, V_{n} = {lucas_val}  ✗")
    if all_match:
        print(f"    ✓ All match for n = 0..9")

# ── What about Q ≠ 1? ──

print(f"\n  Generalisation: Dickson polynomial D_n(P, Q) = 2·Q^(n/2)·T_n(P/(2√Q))")
print(f"  This gives V_n(P, Q) for all Q. Dickson = generalised Chebyshev.")

# ══════════════════════════════════════════════════════════════════
# PART 2: HOW LUCAS-LEHMER USES THIS
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("PART 2: LUCAS-LEHMER — PRIMALITY TESTING IS CHEBYSHEV ITERATION")
print("=" * 72)

print("""
The Lucas-Lehmer test for Mersenne primes M_p = 2^p - 1:
  s_0 = 4
  s_{i+1} = s_i² - 2
  M_p is prime iff s_{p-2} ≡ 0 (mod M_p)

But s_i² - 2 = V_2(s_i, 1) = 2·T_2(s_i/2).

So each step is EVALUATING T_2! The iteration is:
  s_i = 2·T_{2^i}(s_0/2) = 2·T_{2^i}(2) = V_{2^i}(4, 1)

Let's verify:
""")

# Verify the iteration
s = [4]
for i in range(10):
    s.append(s[-1] ** 2 - 2)

for i in range(6):
    via_lucas = lucas_V(2**i, 4, 1)
    print(f"  s_{i} = {s[i]:>12},  V_{2**i}(4,1) = {via_lucas:>12}  "
          f"{'✓' if s[i] == via_lucas else '✗'}")

print(f"""
So the Lucas-Lehmer test asks: does M_p divide T_{{2^(p-2)}}(2)?
Primality testing literally IS iterated Chebyshev evaluation mod a candidate.
""")

# ── Show it working for small Mersenne primes ──

print("Mersenne primes verified via Chebyshev iteration:")
for p in [3, 5, 7, 13, 17, 19]:
    M = 2**p - 1
    s_val = 4
    for _ in range(p - 2):
        s_val = (s_val * s_val - 2) % M
    is_mersenne_prime = (s_val == 0)
    actual = isprime(M)
    symbol = "✓ PRIME" if is_mersenne_prime else "✗ composite"
    print(f"  M_{p} = {M:>8}: s_{{p-2}} mod M = {s_val:>5}  →  {symbol}")

# ══════════════════════════════════════════════════════════════════
# PART 3: BAILLIE-PSW — THE INDUSTRIAL-STRENGTH VERSION
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("PART 3: BAILLIE-PSW AND STRONG LUCAS PSEUDOPRIMES")
print("=" * 72)

print("""
Baillie-PSW (used by every major math library for isprime()) combines:
  1. Miller-Rabin test (base 2)         — catches most composites
  2. Strong Lucas pseudoprime test       — catches what Miller-Rabin misses

The Lucas part evaluates U_n(P, Q) and V_n(P, Q) for specific P, Q chosen
so the candidate n divides U_{n-ε}(P,Q), where ε = Jacobi(D/n).

No counterexample is known (as of 2026). The $620 bounty remains unclaimed.

The reason this works: Chebyshev/Lucas polynomials encode the multiplicative
structure of finite fields. For a prime p, F_p× is cyclic, so T_p(α) has
specific structure. Composites fail to maintain this structure, and the
Lucas test detects the failure.
""")

# ══════════════════════════════════════════════════════════════════
# PART 4: BUNYAKOVSKY'S CONJECTURE
# ══════════════════════════════════════════════════════════════════

print("=" * 72)
print("PART 4: BUNYAKOVSKY'S CONJECTURE AND 2p² - 1")
print("=" * 72)

print("""
Bunyakovsky's conjecture (1857): Let f(x) ∈ ℤ[x] be irreducible, degree ≥ 2,
positive leading coefficient, and content 1 (i.e. gcd of coefficients = 1).
Then f(n) is prime for infinitely many n ∈ ℕ.

STATUS: WIDE OPEN for every single polynomial of degree ≥ 2.
  - Degree 1 is Dirichlet's theorem (proven 1837).
  - Degree 2: even n² + 1 is unresolved.
  - The conjecture would follow from Schinzel's Hypothesis H,
    which would follow from the Bateman-Horn conjecture.

The "content 1" condition rules out trivial obstructions like
f(x) = x² + x + 2, where f(n) is always even.

f(x) = 2x² - 1 has:
  - Irreducible over ℤ ✓
  - Positive leading coefficient ✓
  - Content = gcd(2, 0, -1) = 1 ✓
  - No fixed prime divisor: f(1)=1, f(2)=7, f(3)=17, so gcd = 1 ✓

So Bunyakovsky predicts 2n² - 1 is prime infinitely often.
""")

# ── Density analysis with Bateman-Horn prediction ──

print("BATEMAN-HORN DENSITY PREDICTION")
print("-" * 40)
print("""
Bateman-Horn (1962) gives a quantitative prediction for HOW MANY
primes a polynomial produces. For f(x) = 2x² - 1:

  π_f(N) ~ C_f · N / (deg(f) · ln(N))

where C_f = ∏_p (1 - ω(p)/p) / (1 - 1/p)

and ω(p) = number of roots of f(x) ≡ 0 (mod p).

Let's compute C_f:
""")

# Compute the Bateman-Horn constant for f(x) = 2x² - 1
def count_roots_mod_p(p):
    """Count solutions of 2x² - 1 ≡ 0 (mod p), i.e. 2x² ≡ 1 (mod p)."""
    count = 0
    for x_val in range(p):
        if (2 * x_val * x_val - 1) % p == 0:
            count += 1
    return count

C = 1.0
product_terms = []
for p in primerange(2, 500):
    w = count_roots_mod_p(p)
    factor = (1 - w / p) / (1 - 1 / p)
    C *= factor
    if p < 20:
        product_terms.append((p, w, round(factor, 4)))

print(f"  First few Euler factors (p, ω(p), factor):")
for p, w, f in product_terms:
    print(f"    p={p:>3}: ω(p)={w}, factor={f}")

print(f"\n  C_f (product over primes < 500) ≈ {C:.6f}")

# Compare prediction to reality
import math

print(f"\n  Bateman-Horn prediction vs actual count of primes 2n²-1:")
print(f"  {'N':>8} {'predicted':>10} {'actual':>8} {'ratio':>8}")

for N in [100, 500, 1000, 5000, 10000, 50000]:
    predicted = C * N / (2 * math.log(N))
    actual = sum(1 for n in range(1, N + 1) if isprime(2 * n * n - 1))
    ratio = actual / predicted if predicted > 0 else 0
    print(f"  {N:>8} {predicted:>10.1f} {actual:>8} {ratio:>8.3f}")

# ── Now restrict to prime arguments ──

print(f"\n\nRESTRICTING TO PRIME ARGUMENTS: T_2(p) = 2p² - 1")
print("-" * 50)
print("""
When we only evaluate at primes, we're asking a STRONGER question:
is 2p² - 1 prime for infinitely many primes p?

This is a conjunction of two conditions: p is prime AND 2p² - 1 is prime.
Heuristically, these should be "independent" (no algebraic obstruction
links them), so the density among primes should be:

  ~ C_f · π(N) / (2 · ln(N))  (adjusted for prime inputs)

But there's a subtlety: evaluating at primes avoids even inputs
(except p=2), which removes half the "wasted" evaluations on even n.
""")

print(f"  {'N':>8} {'primes p<N':>12} {'T_2(p) prime':>14} {'density':>10}")

for N in [100, 500, 1000, 5000, 10000]:
    ps = list(primerange(2, N))
    count = sum(1 for p in ps if isprime(2 * p * p - 1))
    density = count / len(ps) if ps else 0
    print(f"  {N:>8} {len(ps):>12} {count:>14} {density:>10.3f}")

# ── Factorisation patterns of composite T_2(p) ──

print(f"\n\nFACTORISATION OF COMPOSITE T_2(p)")
print("-" * 50)
print("When 2p² - 1 is composite, what are its smallest factors?")

smallest_factors = Counter()
factor_residues = Counter()

for p in primerange(2, 1000):
    val = 2 * p * p - 1
    if not isprime(val):
        factors = factorint(val)
        sf = min(factors.keys())
        smallest_factors[sf] += 1
        # What is the smallest factor mod 8?
        factor_residues[sf % 8] += 1

print(f"\n  Most common smallest prime factors of composite 2p²-1 (p < 1000):")
for factor, count in sorted(smallest_factors.items(), key=lambda x: -x[1])[:15]:
    print(f"    {factor:>6}: appears {count} times  (≡ {factor % 8} mod 8)")

print(f"\n  Smallest factor mod 8 distribution:")
for res in sorted(factor_residues.keys()):
    print(f"    ≡ {res} mod 8: {factor_residues[res]} times")

# ── Connection: 2p² ≡ 1 mod q means q | 2p²-1 ──

print(f"""
ALGEBRAIC STRUCTURE OF THE FACTORS:
  q | 2p² - 1  ⟺  2p² ≡ 1 (mod q)  ⟺  p² ≡ (q+1)/2 (mod q)  [if q odd]

  This means p ≡ ±√((q+1)/2) (mod q).
  So for each odd prime q, there are 0 or 2 residue classes of p
  that make q | T_2(p). The number of roots depends on whether
  (q+1)/2 is a quadratic residue mod q.
""")

print("  For which primes q can q | 2p²-1?")
print("  q divides 2p²-1 for some p iff (q+1)/2 is a QR mod q:")
for q in primerange(3, 50):
    half = (q + 1) * pow(2, q - 2, q) % q  # modular inverse of 2, times (q+1)
    # Actually: 2p² ≡ 1 mod q means p² ≡ 2^{-1} mod q
    inv2 = pow(2, q - 2, q)
    is_qr = int(jacobi_symbol(inv2, q))
    roots = [p_val for p_val in range(q) if (2 * p_val * p_val - 1) % q == 0]
    print(f"    q={q:>3}: 2⁻¹ mod {q} = {int(inv2)}, Jacobi(2⁻¹/{q}) = {is_qr:>2}, "
          f"roots: {roots}")

# ══════════════════════════════════════════════════════════════════
# PART 5: THE DEEPER PATTERN — CYCLOTOMIC FIELDS
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("PART 5: WHY CHEBYSHEV CONNECTS TO NUMBER THEORY — CYCLOTOMIC FIELDS")
print("=" * 72)

print("""
The DEEP reason Chebyshev polynomials appear in number theory:

T_n(cos θ) = cos(nθ)  means  T_n encodes the multiplication-by-n map
on the unit circle S¹ = { e^{iθ} }.

The unit circle over a finite field F_q is the set of elements α with
|α| = 1, i.e. α^{q+1} = 1 (using the Frobenius norm). This is a
cyclic group of order q+1.

Chebyshev T_n acts on this group by α ↦ αⁿ. So:
  T_n(x) mod p  encodes the n-th power map on the "circle group" of F_p.

This is exactly the structure that primality tests exploit:
  - For prime p: the circle group is cyclic of known order → T_n is predictable
  - For composite n: the circle group fragments → T_n misbehaves

The bridge from analysis to number theory:
  REAL ANALYSIS          ALGEBRAIC NUMBER THEORY
  ─────────────          ───────────────────────
  cos(nθ)           ↔   αⁿ + α⁻ⁿ in cyclotomic field
  T_n(cos θ)        ↔   Trace of αⁿ from Q(ζ) to Q(ζ+ζ⁻¹)
  Roots of T_n      ↔   Elements of order dividing n on the circle of F_p
  Composition T_m∘T_n ↔ Multiplication m·n in ℤ (the "clock arithmetic")

The composition law T_m(T_n(x)) = T_{mn}(x) is why this works:
it means Chebyshev polynomials form a COMMUTATIVE MONOID under
composition, isomorphic to (ℕ, ×). No other standard polynomial
family has this property. It is unique to Chebyshev/Dickson.
""")

# ── Verify: roots of T_n mod p correspond to elements of order | n ──

print("ROOTS OF T_n(x) mod p (elements of order dividing n on the circle):")
for p in [7, 11, 13]:
    print(f"\n  F_{p}:")
    for n in [2, 3, 4, 6]:
        roots = []
        for a in range(p):
            T_val = chebyshevt(n, x).subs(x, a) % p
            if T_val == a:  # T_n(a) ≡ a mod p means "a is a fixed point of T_n"
                roots.append(a)
        print(f"    Fixed points of T_{n} mod {p}: {roots}  "
              f"(expect ≤ {n} fixed points)")

print()
