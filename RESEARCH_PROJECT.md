# Research Project: Twin Prime Clustering & One-ish Numbers

## Status: Early exploration (Aug 2026)

## The Axiom

We take the **Twin Prime Clustering Conjecture** as an axiom:

**Conjecture (Twin Prime Clustering).** For every positive integer k, there exist infinitely many k-tuples of integers (n_1, n_2, ..., n_k) with n_1 < n_2 < ... < n_k such that for each i in {1, ..., k}, both n_i and n_i + 2 are prime.

**Quantitative form.** Let A = {(a_1, a_1+2), ..., (a_k, a_k+2)} be k twin-pair offsets such that the 2k-element set is admissible. Then the count of n <= N making all 2k values simultaneously prime satisfies:

  pi_A(N) ~ C_A * N / (log N)^(2k)

where C_A = prod over primes q of [(1 - 1/q)^(-2k) * (1 - nu_A(q)/q)], and nu_A(q) is the number of distinct residues mod q hit by the 2k offsets.

### Why it's admissible for all k

Twin pairs at positions {6m, 6m+2} permanently miss residue 1 mod 3, so mod 3 is never an obstruction regardless of k. For primes q >= 5, irregular spacing between pairs avoids covering all residues mod q (greedy/CRT argument).

## One-ish Numbers

**Definition.** A number N is **one-ish of order k** (a **k-one**) if, starting from N+1, there exist k consecutive twin prime pairs. That is, the next k twin prime pairs after N are all present without gaps (Version B — consecutive).

The name comes from the analogy with 1 itself: N is not part of the twin-prime structure, but the structure begins immediately after it. Just as 1 is not prime but the primes begin at 1+1 = 2.

### The hierarchy (filtration)

Under the axiom, we get a nested filtration of the naturals:

  ... c {(k+1)-ones} c {k-ones} c ... c {1-ones} c N

Each level is infinite but exponentially sparser — density ~1/(log x)^(2k). At the limit: **infinity-ones**, numbers from which an infinite sequence of consecutive twin prime pairs begins.

### Compatibility with RH

Twin prime clustering does NOT contradict the Riemann Hypothesis. A run of k consecutive twin prime pairs near x perturbs pi(x) by O(k), while the RH error bound is O(sqrt(x) log x). Since k grows at most like a power of log x (by density predictions), the two live at completely different scales.

RH controls global regularity; twin prime clustering is local correlation. They are allies, not adversaries — GRH is in fact used as input to the sieve methods that approach bounded prime gaps.

## Open threads

### 1. Self-similar structure of the hierarchy
Do k-ones themselves cluster? If k-ones behave pseudo-randomly at density 1/(log x)^(2k), then Bateman-Horn heuristics predict "twin k-ones" at density 1/(log x)^(4k) — a hierarchy of hierarchies. Every conjecture about primes (twin primes, Goldbach, prime number theorem) may have a shadow at each level.

### 2. The infinity-ones
The most mysterious objects. Their existence is guaranteed by the axiom but their distribution is opaque. These are the "fixed points" of the hierarchy — numbers where twin-prime structure never runs out. Are these the right objects to study for potential tension with RH?

### 3. Category-theoretic flavour
A k-one is characterised not by what it IS but by what unfolds after it — its significance is entirely relational, forward-looking. This echoes Yoneda: an object characterised by its morphisms out. Is there a formal category here? Objects = naturals, morphisms = twin-prime witnessing chains?

### 4. Connection to existing project
The constant C_A for a cluster configuration is computed by the same Bateman-Horn machinery as in `galois_beat_euler.py` and `prime_friendliness.py`. The "shield" concept (which small primes are blocked) determines which cluster shapes are most one-ish-friendly.
