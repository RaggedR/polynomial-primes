import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from sympy import isprime, primerange

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'

primes_list = list(primerange(2, 10000))

polys = [
    (r'$6p^2 - 6p + 1$  (Bernoulli $B_2$)',   lambda p: 6*p*p - 6*p + 1,   '#e63946'),
    (r'$2p^2 - 1$  (Chebyshev $T_2$)',         lambda p: 2*p*p - 1,         '#457b9d'),
    (r'$p^2 - p - 1$  (golden ratio)',          lambda p: p*p - p - 1,       '#2a9d8f'),
    (r'$8p^2 - 31$  (Meixner–Pollaczek $P_7$)',lambda p: 8*p*p - 31,        '#e9c46a'),
    (r'$4p^2 - 3$  (Chebyshev $T_3$ factor)',   lambda p: 4*p*p - 3,        '#9b5de5'),
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
    ax.plot(ps, sm, color=color, linewidth=2.2, label=label)
    print(f"{label[:40]:>40}  {total}/{len(primes_list)} = {total/len(primes_list):.3f}")

ax.set_xlabel('Prime $p$', fontsize=13)
ax.set_ylabel('Prime density', fontsize=13)
ax.set_title('Top 5 quadratic prime generators at prime arguments',
             fontsize=15, fontweight='bold')
ax.legend(fontsize=10, loc='upper right')
ax.grid(True, alpha=0.3)
ax.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig('/Users/robin/git/primes/top5.png', dpi=200, bbox_inches='tight')
plt.savefig('/Users/robin/git/primes/top5.pdf', bbox_inches='tight')
print("Done")
