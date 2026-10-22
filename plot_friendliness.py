"""
Plot prime-friendliness of Chebyshev polynomials T_{2^k}.
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
from sympy import chebyshevt, Symbol, Poly, QQ, isprime, primerange

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'

x = Symbol('x')

# ── Precompute data ──────────────────────────────────────────────

sieve_primes = list(primerange(2, 200))

def count_roots_mod_p(f_poly, p):
    coeffs = Poly(f_poly, x, domain=QQ).all_coeffs()
    count = 0
    for x_val in range(p):
        val = 0
        for c in coeffs:
            val = (val * x_val + int(c)) % p
        if val == 0:
            count += 1
    return count

def bateman_horn_constant(f_poly):
    C = 1.0
    for p in sieve_primes:
        w = count_roots_mod_p(f_poly, p)
        C *= (1 - w / p) / (1 - 1 / p)
    return C

def eval_T(poly_expr, val):
    return int(poly_expr.subs(x, val))

# Compute for T_1, T_2, T_4, T_8, T_16
ks = [0, 1, 2, 3, 4]
labels = [r'$T_1 = x$', r'$T_2 = 2x^2 - 1$', r'$T_4$', r'$T_8$', r'$T_{16}$']
colors = ['#888888', '#e63946', '#457b9d', '#2a9d8f', '#e9c46a']

polys = {2**k: chebyshevt(2**k, x) for k in ks}
constants = {2**k: bateman_horn_constant(polys[2**k]) for k in ks}

N_max = 1000
Ns = list(range(2, N_max + 1))

# Count primes cumulatively
print("Computing prime counts...")
cumulative = {}
for k in ks:
    n = 2**k
    T = polys[n]
    counts = []
    total = 0
    for m in Ns:
        val = eval_T(T, m)
        if abs(val) > 1 and isprime(abs(val)):
            total += 1
        counts.append(total)
    cumulative[n] = counts
    print(f"  T_{n}: {total} primes up to {N_max}")

# ── Figure 1: Cumulative prime counts vs Bateman-Horn prediction ──

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Top-left: cumulative prime counts
ax = axes[0, 0]
for k, label, color in zip(ks, labels, colors):
    n = 2**k
    ax.plot(Ns, cumulative[n], color=color, label=label, linewidth=1.8)
    # Bateman-Horn prediction
    deg = n
    C = constants[n]
    predicted = [C * m / (deg * math.log(m)) for m in Ns]
    ax.plot(Ns, predicted, color=color, linestyle='--', alpha=0.5, linewidth=1)

ax.set_xlabel('$N$')
ax.set_ylabel('Primes found')
ax.set_title('Cumulative prime count (solid) vs Bateman–Horn (dashed)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Top-right: running density (primes / N)
ax = axes[0, 1]
for k, label, color in zip(ks, labels, colors):
    n = 2**k
    density = [cumulative[n][i] / (i + 2) for i in range(len(Ns))]
    # Smooth with moving average
    window = 50
    smoothed = np.convolve(density, np.ones(window)/window, mode='valid')
    ax.plot(Ns[window-1:], smoothed, color=color, label=label, linewidth=1.8)

ax.set_xlabel('$N$')
ax.set_ylabel('Density (primes found / $N$)')
ax.set_title('Running prime density (smoothed)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Bottom-left: C_f and C_f/degree bar chart
ax = axes[1, 0]
bar_x = np.arange(len(ks))
bar_width = 0.35

cf_vals = [constants[2**k] for k in ks]
cf_over_deg = [constants[2**k] / (2**k) for k in ks]

bars1 = ax.bar(bar_x - bar_width/2, cf_vals, bar_width,
               color=[c for c in colors], alpha=0.7, label=r'$C_f$ (shield strength)')
bars2 = ax.bar(bar_x + bar_width/2, cf_over_deg, bar_width,
               color=[c for c in colors], alpha=1.0, edgecolor='black',
               linewidth=1.2, label=r'$C_f / \deg$ (effective rate)')

ax.set_xticks(bar_x)
ax.set_xticklabels([f'$T_{{{2**k}}}$' for k in ks])
ax.set_ylabel('Value')
ax.set_title(r'Shield strength $C_f$ vs effective rate $C_f / \deg$')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis='y')

# Annotate the bars
for bar, val in zip(bars1, cf_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
            f'{val:.2f}', ha='center', va='bottom', fontsize=8)
for bar, val in zip(bars2, cf_over_deg):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
            f'{val:.2f}', ha='center', va='bottom', fontsize=8)

# Bottom-right: shield visualisation (blocked primes heatmap)
ax = axes[1, 1]

small_primes = list(primerange(2, 60))
shield = np.zeros((len(small_primes), len(ks)))

for j, k in enumerate(ks):
    n = 2**k
    T = polys[n]
    for i, q in enumerate(small_primes):
        w = count_roots_mod_p(T, q)
        shield[i, j] = w

# Custom colormap: 0 = blocked (dark blue), >0 = exposed (red scale)
from matplotlib.colors import ListedColormap, BoundaryNorm
cmap = ListedColormap(['#1d3557', '#f4a261', '#e76f51', '#d62828', '#6a040f'])
bounds = [-0.5, 0.5, 1.5, 2.5, 4.5, 16.5]
norm = BoundaryNorm(bounds, cmap.N)

im = ax.imshow(shield, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')
ax.set_xticks(range(len(ks)))
ax.set_xticklabels([f'$T_{{{2**k}}}$' for k in ks])
ax.set_yticks(range(0, len(small_primes), 2))
ax.set_yticklabels([str(small_primes[i]) for i in range(0, len(small_primes), 2)])
ax.set_ylabel('Prime $q$')
ax.set_title('Small prime shield (dark = blocked)')
cb = fig.colorbar(im, ax=ax, ticks=[0, 1, 2, 4, 8, 16], shrink=0.8)
cb.set_label(r'$\omega(q)$ = roots mod $q$')

fig.suptitle('Prime-Friendliness of Chebyshev Polynomials $T_{2^k}$',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('/Users/robin/git/primes/prime_friendliness.png', dpi=200, bbox_inches='tight')
plt.savefig('/Users/robin/git/primes/prime_friendliness.pdf', bbox_inches='tight')
print("\nSaved to prime_friendliness.png and .pdf")
