import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from sympy import isprime

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'

N = 5000

counts_f1 = []
counts_f2 = []
t1 = t2 = 0

for n in range(2, N + 1):
    v1 = 2 * n * n - 1
    v2 = 8 * n * n - 31
    if v1 > 1 and isprime(v1):
        t1 += 1
    if v2 > 1 and isprime(abs(v2)):
        t2 += 1
    counts_f1.append(t1)
    counts_f2.append(t2)

ns = list(range(2, N + 1))

fig, ax = plt.subplots(figsize=(10, 6))

window = 100
density_f1 = [counts_f1[i] / (i + 2) for i in range(len(ns))]
density_f2 = [counts_f2[i] / (i + 2) for i in range(len(ns))]
sm_f1 = np.convolve(density_f1, np.ones(window)/window, mode='valid')
sm_f2 = np.convolve(density_f2, np.ones(window)/window, mode='valid')

ax.plot(ns[window-1:], sm_f1, color='#e63946', linewidth=2.2,
        label=r'$2n^2 - 1$  (Chebyshev $T_2$)')
ax.plot(ns[window-1:], sm_f2, color='#457b9d', linewidth=2.2,
        label=r'$8n^2 - 31$  (Meixner–Pollaczek $P_7$ factor)')

ax.set_xlabel('$n$', fontsize=13)
ax.set_ylabel('Prime density', fontsize=13)
ax.set_title('Running prime density of two quadratics', fontsize=15, fontweight='bold')
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
ax.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig('/Users/robin/git/primes/density_big.png', dpi=200, bbox_inches='tight')
plt.savefig('/Users/robin/git/primes/density_big.pdf', bbox_inches='tight')
print("Done")
