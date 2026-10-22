"""
Constructing quintics to beat Euler using class field theory.

Strategy 1: CYCLIC QUINTICS (G = ℤ/5)
  - Splitting controlled by conductor: q splits completely iff q ≡ 1 mod f
  - 80% of primes blocked (ω=0), but 20% have ω=5 (harsh penalty)
  - Conductor 11: first split prime is q=23 → all primes < 23 blocked

Strategy 2: DIHEDRAL/FROBENIUS QUINTICS
  - F_20 group: gentler ω=1 breakthroughs, but only 20% blocked
  - Need the blocked 20% to be the SMALL primes

Strategy 3: QUINTICS FROM NUMBER FIELDS WITH SMALL DISCRIMINANT
  - Analogue of Heegner: find quintic fields with special properties

We use known quintic polynomials generating specific number fields.
"""

from sympy import Symbol, Poly, QQ, factor, isprime, primerange, factorint, cyclotomic_poly
from math import gcd, log
from functools import reduce

x = Symbol('x')
sieve_primes = list(primerange(2, 500))
eval_primes = list(primerange(2, 10000))

def count_roots_mod_p(coeffs, p):
    count = 0
    for xv in range(p):
        val = 0
        for c in coeffs:
            val = (val * xv + c) % p
        if val == 0:
            count += 1
    return count

def bateman_horn_fast(coeffs):
    C = 1.0
    for p in sieve_primes:
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

def full_analysis(name, coeffs, deg=5):
    C = bateman_horn_fast(coeffs)
    Cd = C / deg

    # Blocked primes
    blocked = []
    omega_dist = {}
    for q in primerange(2, 80):
        w = count_roots_mod_p(coeffs, q)
        if w == 0:
            blocked.append(q)
        omega_dist[q] = w

    # Prime density
    count = 0
    for p in eval_primes:
        val = abs(eval_poly(coeffs, p))
        if val > 1 and isprime(val):
            count += 1
    density = count / len(eval_primes)

    # Irreducibility
    sym_poly = sum(c * x**i for i, c in enumerate(reversed(coeffs)))
    fac = factor(sym_poly, domain=QQ)
    is_irr = str(fac) == str(sym_poly)

    return {
        'name': name, 'coeffs': coeffs, 'C': C, 'Cd': Cd,
        'blocked': blocked, 'omega': omega_dist,
        'count': count, 'density': density, 'irr': is_irr,
        'poly': sym_poly
    }

# ══════════════════════════════════════════════════════════════════
# EULER BASELINE
# ══════════════════════════════════════════════════════════════════

euler = full_analysis("Euler p²+p+41", [1, 1, 41], deg=2)
print("=" * 75)
print(f"EULER BASELINE: C_f = {euler['C']:.3f}, C_f/deg = {euler['Cd']:.4f}, "
      f"density = {euler['density']:.3f}")
print(f"  Blocked: {euler['blocked']}")
print(f"  TARGET: C_f > {euler['Cd'] * 5:.2f} for a quintic")
print("=" * 75)

# ══════════════════════════════════════════════════════════════════
# STRATEGY 1: CYCLIC QUINTICS FROM CYCLOTOMIC FIELDS
# ══════════════════════════════════════════════════════════════════

print(f"\n{'='*75}")
print("STRATEGY 1: CYCLIC QUINTICS (G = ℤ/5)")
print("="*75)
print("""
For conductor f, a cyclic quintic is a factor of Φ_f(x) when 5 | φ(f).
Prime q has ω(q) = 5 if q ≡ 1 mod f, else ω(q) = 0 (ignoring ramification).

The smallest conductors with cyclic quintic subfields:
  f = 11: φ(11) = 10, contains ℤ/5. First split prime: q = 23.
  f = 25: φ(25) = 20, contains ℤ/5. First split prime: q = 101.
  f = 31: φ(31) = 30, contains ℤ/5. First split prime: q = ?
""")

# The cyclic quintic subfield of ℚ(ζ_11) has a known minimal polynomial.
# Gauss showed: x⁵ + x⁴ - 4x³ - 3x² + 3x + 1 generates the maximal real
# subfield of ℚ(ζ_11), which is cyclic of degree 5.

cyclic_quintics = [
    # (name, coeffs [highest to lowest])
    ("ℚ(ζ₁₁)⁺ (conductor 11)",
     [1, 1, -4, -3, 3, 1]),

    # Quintic subfield of ℚ(ζ_25)
    ("ℚ(ζ₂₅)⁺ (conductor 25)",
     [1, 0, -5, 0, 5, -1]),  # Φ_25 has a quintic factor

    # Quintic subfield of ℚ(ζ_31)
    ("ℚ(ζ₃₁) quintic (conductor 31)",
     [1, 1, -12, -21, 1, 5]),  # known quintic subfield

    # Emma Lehmer's quintic (conductor 25, related)
    ("Lehmer quintic (a=1)",
     [1, 1, -4, -3, 3, 1]),  # Same as conductor 11 case

    # Conductor 41
    ("ℚ(ζ₄₁) quintic (conductor 41)",
     [1, 1, -16, -5, 21, -1]),  # quintic subfield

    # Simple cyclic: x⁵ - 5x³ + 5x - 1 = ... let me use known ones
    ("x⁵-2 (G = F_20, non-abelian!)",
     [1, 0, 0, 0, 0, -2]),

    ("x⁵-x-1 (G = S_5)",
     [1, 0, 0, 0, -1, -1]),

    ("x⁵+x+1 (has x²+x+1 factor)",
     [1, 0, 0, 0, 1, 1]),

    ("x⁵-5x+12",
     [1, 0, 0, 0, -5, 12]),

    ("x⁵-x²+1",
     [1, 0, 0, -1, 0, 1]),
]

results = []
for name, coeffs in cyclic_quintics:
    r = full_analysis(name, coeffs)
    results.append(r)

    marker = " ★ BEATS EULER" if r['Cd'] > euler['Cd'] else ""
    irr = "IRR" if r['irr'] else "RED"
    print(f"\n  {name}")
    print(f"    {r['poly']}")
    print(f"    {irr}  C_f={r['C']:.3f}  C_f/5={r['Cd']:.4f}  "
          f"density={r['density']:.3f}  ({r['count']}/{len(eval_primes)}){marker}")
    print(f"    Blocked < 80: {r['blocked']}")
    # Show ω distribution for small primes
    omega_str = " ".join(f"{q}:{r['omega'][q]}" for q in list(primerange(2, 30)))
    print(f"    ω(q): {omega_str}")


# ══════════════════════════════════════════════════════════════════
# STRATEGY 2: CONSTRUCT FROM CONDUCTOR THEORY
# ══════════════════════════════════════════════════════════════════

print(f"\n{'='*75}")
print("STRATEGY 2: SEARCH NEAR KNOWN CYCLIC QUINTICS")
print("="*75)
print("Perturb the conductor-11 quintic to find better variants\n")

base = [1, 1, -4, -3, 3, 1]  # conductor 11 quintic

best = []
for da in range(-3, 4):
    for db in range(-8, 9):
        for dc in range(-12, 13):
            for dd in range(-12, 13):
                for de in range(-20, 21):
                    coeffs = [1, base[1]+da, base[2]+db, base[3]+dc,
                              base[4]+dd, base[5]+de]
                    # Quick screen
                    blocked = sum(1 for p in [2,3,5,7,11,13,17,19,23]
                                  if count_roots_mod_p(coeffs, p) == 0)
                    if blocked < 7:
                        continue
                    C = bateman_horn_fast(coeffs)
                    if C / 5 > euler['Cd'] * 0.9:
                        sym = sum(c * x**i for i, c in enumerate(reversed(coeffs)))
                        fac = factor(sym, domain=QQ)
                        is_irr = str(fac) == str(sym)
                        if is_irr:
                            best.append((C, C/5, coeffs, sym))

best.sort(key=lambda x: -x[0])

print(f"  Found {len(best)} irreducible quintics with C_f/5 > {euler['Cd']*0.9:.2f}")
print(f"\n  {'C_f':>7} {'C/5':>6} {'Polynomial'}")
print(f"  {'-'*60}")
for C, Cd, coeffs, sym in best[:20]:
    marker = " ★" if Cd > euler['Cd'] else ""
    print(f"  {C:>7.3f} {Cd:>6.3f} {sym}{marker}")


# ══════════════════════════════════════════════════════════════════
# EVALUATE BEST CANDIDATES
# ══════════════════════════════════════════════════════════════════

if best:
    print(f"\n{'='*75}")
    print("BEST CANDIDATES: FULL EVALUATION")
    print("="*75)

    # Take top 5 and compute actual density
    for C, Cd, coeffs, sym in best[:5]:
        r = full_analysis(str(sym), coeffs)
        marker = " ★ BEATS EULER" if r['density'] > euler['density'] else ""
        print(f"\n  {sym}")
        print(f"    C_f={r['C']:.3f}  C_f/5={r['Cd']:.4f}  "
              f"density={r['density']:.3f}  ({r['count']}/{len(eval_primes)}){marker}")
        print(f"    Blocked: {r['blocked']}")
        omega_str = " ".join(f"{q}:{r['omega'][q]}" for q in list(primerange(2, 30)))
        print(f"    ω(q): {omega_str}")

    # Compare everything
    print(f"\n{'='*75}")
    print("FINAL COMPARISON")
    print("="*75)
    print(f"\n  {'Name':>35} {'deg':>4} {'C_f':>7} {'C/d':>6} {'density':>8} {'beat?':>6}")
    print(f"  {'-'*70}")
    print(f"  {'Euler p²+p+41':>35} {2:>4} {euler['C']:>7.3f} {euler['Cd']:>6.3f} "
          f"{euler['density']:>8.3f}")

    for C, Cd, coeffs, sym in best[:10]:
        count = 0
        for p in eval_primes:
            val = abs(eval_poly(coeffs, p))
            if val > 1 and isprime(val):
                count += 1
        dens = count / len(eval_primes)
        beat = "★" if dens > euler['density'] else ""
        short = str(sym)
        if len(short) > 35:
            short = short[:32] + "..."
        print(f"  {short:>35} {5:>4} {C:>7.3f} {Cd:>6.3f} {dens:>8.3f} {beat:>6}")
