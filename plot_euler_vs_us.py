import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from sympy import isprime, primerange
import math

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'

primes_list = list(primerange(2, 10000))

# "Random" baseline: probability a number near p² is prime ≈ 1/(2·ln(p))
# We compute the actual density of primes among odd numbers near p²

polys = [
    (r"$p^2 + p + 41$  (Euler, 1772)",     lambda p: p*p + p + 41,   '#e63946'),
    (r"$6p^2 - 6p + 1$  (Bernoulli $B_2$)",lambda p: 6*p*p - 6*p + 1,'#457b9d'),
    (r"$2p^2 - 1$  (Chebyshev $T_2$)",     lambda p: 2*p*p - 1,      '#2a9d8f'),
]

fig, ax = plt.subplots(figsize=(11, 6.5))
window = 50

for label, f, color in polys:
    counts = []
    total = 0
    for p in primes_list:
        val = f(p)
        if abs(val) > 1 and isprime(abs(val)):
            total += 1
        counts.append(total)

    density = [counts[i] / (i + 1) for i in range(len(primes_list))]
    sm = np.convolve(density, np.ones(window)/window, mode='valid')
    ps = [primes_list[i] for i in range(window - 1, len(primes_list))]
    ax.plot(ps, sm, color=color, linewidth=2.5, label=label)

# Random baseline: 1/ln(p²) = 1/(2·ln(p))
ps_all = [primes_list[i] for i in range(window - 1, len(primes_list))]
random_density = [1 / (2 * math.log(p)) for p in ps_all]
ax.plot(ps_all, random_density, color='#888888', linewidth=2, linestyle='--',
        label=r'Random odd number near $p^2$:  $1/\ln p^2$')

ax.set_xlabel('Prime $p$', fontsize=14)
ax.set_ylabel('Prime density', fontsize=14)
ax.set_title('Euler vs Bernoulli vs Chebyshev vs random',
             fontsize=16, fontweight='bold')
ax.legend(fontsize=11.5, loc='upper right')
ax.grid(True, alpha=0.3)
ax.tick_params(labelsize=12)

plt.tight_layout()
plt.savefig('/Users/robin/git/primes/euler_vs_us.jpg', dpi=200, bbox_inches='tight')
plt.savefig('/Users/robin/git/primes/euler_vs_us.pdf', bbox_inches='tight')
print("Done")
