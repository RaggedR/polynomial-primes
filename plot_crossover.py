import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from sympy import isprime, primerange

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'

primes_list = list(primerange(2, 50001))

euler_f = lambda p: p*p + p + 41
rival_f = lambda p: p*p + p + 21377

# Cumulative prime counts
euler_cum = []
rival_cum = []
e_tot = r_tot = 0
for p in primes_list:
    ev = euler_f(p)
    rv = rival_f(p)
    if ev > 1 and isprime(ev):
        e_tot += 1
    if rv > 1 and isprime(rv):
        r_tot += 1
    euler_cum.append(e_tot)
    rival_cum.append(r_tot)

euler_cum = np.array(euler_cum, dtype=float)
rival_cum = np.array(rival_cum, dtype=float)
p_arr = np.array(primes_list)

# Compute running density (cumulative)
euler_dens = euler_cum / np.arange(1, len(primes_list) + 1)
rival_dens = rival_cum / np.arange(1, len(primes_list) + 1)

# Smooth with window
window = 80
euler_sm = np.convolve(euler_dens, np.ones(window)/window, mode='valid')
rival_sm = np.convolve(rival_dens, np.ones(window)/window, mode='valid')
ps = p_arr[window - 1:]

fig, ax = plt.subplots(figsize=(11, 6.5))

ax.plot(ps, euler_sm, color='#457b9d', linewidth=2.5,
        label=r'$p^2 + p + 41$  (Euler, $\Delta = -163$)')
ax.plot(ps, rival_sm, color='#e63946', linewidth=2.5,
        label=r'$p^2 + p + 21377$  ($\Delta = -85507$)')

# Find crossover point
diff = rival_sm - euler_sm
crossings = np.where(np.diff(np.sign(diff)))[0]
# Take the last major crossing (the definitive one)
if len(crossings) > 0:
    cx = crossings[-1]
    cross_p = int(ps[cx])
    cross_dens = float(euler_sm[cx])
    ax.axvline(cross_p, color='#888888', linewidth=1, linestyle=':', alpha=0.7)
    ax.annotate(f'crossover\n$p \\approx {cross_p:,}$',
                xy=(cross_p, cross_dens), xytext=(cross_p + 3000, cross_dens + 0.015),
                fontsize=12, color='#555555',
                arrowprops=dict(arrowstyle='->', color='#888888', lw=1.2))

# Zoom: set y-limits to make the gap visible
y_min = min(euler_sm.min(), rival_sm.min()) - 0.01
y_max = max(euler_sm[0], rival_sm[0]) + 0.01
ax.set_ylim(y_min, y_max)

ax.set_xlabel('Prime $p$', fontsize=14)
ax.set_ylabel('Prime density (cumulative)', fontsize=14)
ax.set_title(r'The Crossover: $\Delta=-85507$ overtakes Euler',
             fontsize=16, fontweight='bold')
ax.legend(fontsize=12, loc='upper right')
ax.grid(True, alpha=0.3)
ax.tick_params(labelsize=12)

plt.tight_layout()
plt.savefig('/Users/robin/git/primes/crossover.jpg', dpi=200, bbox_inches='tight')
plt.savefig('/Users/robin/git/primes/crossover.pdf', bbox_inches='tight')
print("Done")
