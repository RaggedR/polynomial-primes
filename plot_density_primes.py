import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from sympy import isprime, primerange

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'

primes_list = list(primerange(2, 10000))

counts_f1 = []
counts_f2 = []
t1 = t2 = 0

for p in primes_list:
    v1 = 2 * p * p - 1
    v2 = 8 * p * p - 31
    if v1 > 1 and isprime(v1):
        t1 += 1
    if v2 > 1 and isprime(abs(v2)):
        t2 += 1
    counts_f1.append(t1)
    counts_f2.append(t2)

print(f"Primes p < 10000: {len(primes_list)}")
print(f"  2p²-1  prime: {t1}/{len(primes_list)} = {t1/len(primes_list):.3f}")
print(f"  8p²-31 prime: {t2}/{len(primes_list)} = {t2/len(primes_list):.3f}")

fig, ax = plt.subplots(figsize=(10, 6))

window = 50
density_f1 = [counts_f1[i] / (i + 1) for i in range(len(primes_list))]
density_f2 = [counts_f2[i] / (i + 1) for i in range(len(primes_list))]
sm_f1 = np.convolve(density_f1, np.ones(window)/window, mode='valid')
sm_f2 = np.convolve(density_f2, np.ones(window)/window, mode='valid')
ps = [primes_list[i] for i in range(window - 1, len(primes_list))]

ax.plot(ps, sm_f1, color='#e63946', linewidth=2.2,
        label=r'$2p^2 - 1$  (Chebyshev $T_2$)')
ax.plot(ps, sm_f2, color='#457b9d', linewidth=2.2,
        label=r'$8p^2 - 31$  (Meixner–Pollaczek $P_7$ factor)')

ax.set_xlabel('Prime $p$', fontsize=13)
ax.set_ylabel('Prime density', fontsize=13)
ax.set_title('Running prime density at prime arguments', fontsize=15, fontweight='bold')
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
ax.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig('/Users/robin/git/primes/density_primes.png', dpi=200, bbox_inches='tight')
plt.savefig('/Users/robin/git/primes/density_primes.pdf', bbox_inches='tight')
print("Done")
