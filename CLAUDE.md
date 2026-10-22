# CLAUDE.md — primes project

## What this project is

An exploratory number theory project investigating **prime values of classical polynomial families** — Chebyshev, Hermite, Meixner-Pollaczek, Bernoulli, Euler, cyclotomic, Dickson, Gegenbauer, Legendre, Laguerre, Fibonacci — evaluated at prime arguments.

## Key findings

1. **Chebyshev T₂(p) = 2p²−1** produces primes at high density (~20% at prime arguments p < 10000). This connects to Lucas sequences and primality testing (Lucas-Lehmer is iterated Chebyshev evaluation). T_n is irreducible over ℚ iff n is a power of 2.

2. **Bernoulli B₂ factor 6p²−6p+1** beats Chebyshev (C_f/deg = 2.06 vs 1.84). The Meixner-Pollaczek P₇ factor 8p²−31 also beats it (C_f = 3.88).

3. **Euler's p²+p+41 is NOT the Bateman-Horn optimal quadratic.** Discriminant Δ=−85507 (= −37·2311) gives C_f/2 = 3.48 vs Euler's 3.34. The mechanism: 37 is ramified (neutral), and ALL primes from 2 to 59 are blocked. The crossover is slow — the polynomials trade leads up to p < 50000. Euler's Δ=−163 remains optimal among class-number-1 discriminants (Heegner), but composite discriminants exploit ramification to extend the shield.

4. **No polynomial of degree ≥ 3 can beat quadratics.** The "degree tax" (dividing C_f by deg) is devastating. A cubic needs to block all primes up to ~251 (vs 37-61 for quadratics); quartics need ~1759; quintics ~11789. The best cyclic cubic found (conductor 4851139, blocking all primes up to 131) gets C_f/3 = 2.88 — still 14% behind Euler. The structural obstruction: E[Bateman-Horn factor] = 1 for ANY Galois group (a theorem), so the only lever is which SPECIFIC small primes are blocked.

5. **The prime-friendliness of a polynomial is controlled by:** (a) which small primes are blocked (quadratic reciprocity / Chebotarev), (b) how harsh the breakthrough penalty is (determined by Galois group fixed-point statistics), (c) the degree. Open: is there a theoretical framework (analogous to Heegner) for optimizing composite discriminants?

## Files

- `explore.py` — initial exploration of T_n(p), H_n(p)
- `lucas_bunyakovsky.py` — deep dive: Chebyshev=Lucas identity, Lucas-Lehmer, Bunyakovsky, Bateman-Horn
- `bunyakovsky.tex/.pdf` — expository paper on Bunyakovsky's conjecture
- `irreducibility.py` — proves T_n irreducible iff n = 2^k
- `prime_friendliness.py` — Bateman-Horn constants for T_{2^k}
- `plot_friendliness.py` — 4-panel plot: counts, density, C_f bar chart, shield heatmap
- `hermite_density.py` — Hermite polynomial analysis (H_n(p) is never prime; cores analysed)
- `mp_density.py` — Meixner-Pollaczek analysis; found 8p²−31 beats T₂
- `survey_families.py` — systematic survey of 10 polynomial families; found Bernoulli B₂ wins among orthogonal polys
- `euler41.py` — comparison with Euler's p²+p+41; Euler wins decisively
- `quintic_cft.py` — class field theory approach to beating Euler with quintics; proves it's structurally impossible
- `galois_beat_euler.py` — Galois theory attack: shield bounds by group, exhaustive discriminant search, cyclic cubic conductor search; found Δ=−85507 beats Euler on C_f
- `chebyshev_decomposition.py` — additive decomposition into "Chebyshev primes" (non-unique)
- `generate_large.py` — practical prime generation benchmarks
- Various plot scripts and outputs (PNG, PDF, JPG)

## Open threads

- **What makes composite discriminants optimal?** Δ=−85507=−37·2311 beats Euler by exploiting ramification at 37 to extend the shield. Is there a theoretical framework predicting which composite discriminants maximize C_f? The interplay of ramification and quadratic reciprocity at composite moduli is unexplored.
- **Deeper cubic search**: best conductor found blocks primes to 137 with q ≤ 10^7. Searching q ≤ 10^10 might find conductors blocking to ~200+, but C_f/3 still needs to exceed 3.48 — likely impossible.
- Connection to **modular forms**: the j-invariant / class polynomial angle could give a deeper understanding of why certain discriminants produce better shields.
- The **function field analogue** of Bunyakovsky is a theorem (via Weil). Could the proof strategy inform the number field case?

## Dependencies

- Python 3, sympy, matplotlib, numpy, PIL
- tectonic (for LaTeX compilation)
