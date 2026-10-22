"""
Head-to-head: 8x² - 31 vs 2x² - 1 prime density.
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from sympy import isprime, primerange

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'

N = 5000

# ── Compute cumulative prime counts for all n ────────────────────

def f1(n): return 2 * n * n - 1       # T_2 = Chebyshev
def f2(n): return 8 * n * n - 31      # MP_7 factor = Meixner-Pollaczek

counts_f1 = []
counts_f2 = []
total_f1 = 0
total_f2 = 0

for n in range(2, N + 1):
    v1 = f1(n)
    v2 = f2(n)
    if v1 > 1 and isprime(v1):
        total_f1 += 1
    if v2 > 1 and isprime(abs(v2)):
        total_f2 += 1
    counts_f1.append(total_f1)
    counts_f2.append(total_f2)

ns = list(range(2, N + 1))

print(f"Total primes up to n={N}:")
print(f"  2n²-1:  {total_f1}")
print(f"  8n²-31: {total_f2}")

# ── Also at prime arguments only ─────────────────────────────────

primes_list = list(primerange(2, N))
prime_counts_f1 = []
prime_counts_f2 = []
t1 = t2 = 0
prime_ns = []

for p in primes_list:
    v1 = f1(p)
    v2 = f2(p)
    if v1 > 1 and isprime(v1):
        t1 += 1
    if v2 > 1 and isprime(abs(v2)):
        t2 += 1
    prime_counts_f1.append(t1)
    prime_counts_f2.append(t2)
    prime_ns.append(p)

print(f"\nAt prime arguments p < {N}:")
print(f"  2p²-1:  {t1}/{len(primes_list)} = {t1/len(primes_list):.3f}")
print(f"  8p²-31: {t2}/{len(primes_list)} = {t2/len(primes_list):.3f}")

# ── Plot ─────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

red = '#e63946'
blue = '#457b9d'

# Top-left: cumulative prime counts (all n)
ax = axes[0, 0]
ax.plot(ns, counts_f1, color=red, linewidth=1.8, label=r'$2n^2 - 1$ (Chebyshev $T_2$)')
ax.plot(ns, counts_f2, color=blue, linewidth=1.8, label=r'$8n^2 - 31$ (Meixner–Pollaczek)')
ax.set_xlabel('$n$')
ax.set_ylabel('Primes found')
ax.set_title('Cumulative prime count (all $n$)')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Top-right: running density (all n)
ax = axes[0, 1]
window = 100
density_f1 = [counts_f1[i] / (i + 2) for i in range(len(ns))]
density_f2 = [counts_f2[i] / (i + 2) for i in range(len(ns))]
sm_f1 = np.convolve(density_f1, np.ones(window)/window, mode='valid')
sm_f2 = np.convolve(density_f2, np.ones(window)/window, mode='valid')
ax.plot(ns[window-1:], sm_f1, color=red, linewidth=1.8, label=r'$2n^2 - 1$')
ax.plot(ns[window-1:], sm_f2, color=blue, linewidth=1.8, label=r'$8n^2 - 31$')
ax.set_xlabel('$n$')
ax.set_ylabel('Density')
ax.set_title('Running prime density (smoothed, all $n$)')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Bottom-left: cumulative at prime arguments
ax = axes[1, 0]
ax.plot(prime_ns, prime_counts_f1, color=red, linewidth=1.8, label=r'$2p^2 - 1$')
ax.plot(prime_ns, prime_counts_f2, color=blue, linewidth=1.8, label=r'$8p^2 - 31$')
ax.set_xlabel('Prime $p$')
ax.set_ylabel('Primes found')
ax.set_title('Cumulative prime count (prime arguments only)')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Bottom-right: ratio f2/f1 over time
ax = axes[1, 1]
ratio = [counts_f2[i] / counts_f1[i] if counts_f1[i] > 0 else 0
         for i in range(len(ns))]
sm_ratio = np.convolve(ratio, np.ones(window)/window, mode='valid')
ax.plot(ns[window-1:], sm_ratio, color='#2a9d8f', linewidth=1.8)
ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
ax.axhline(y=3.876/3.675, color='orange', linestyle='--', alpha=0.7,
           label=f'$C_{{8x^2-31}} / C_{{2x^2-1}}$ = {3.876/3.675:.3f}')
ax.set_xlabel('$n$')
ax.set_ylabel('Ratio')
ax.set_title(r'Ratio: primes from $8n^2\!-31$ / primes from $2n^2\!-1$')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_ylim(0.9, 1.2)

fig.suptitle(r'Chebyshev $T_2$ vs Meixner–Pollaczek: which quadratic produces more primes?',
             fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('/Users/robin/git/primes/duel.png', dpi=200, bbox_inches='tight')
plt.savefig('/Users/robin/git/primes/duel.pdf', bbox_inches='tight')
print("\nSaved to duel.png and duel.pdf")
