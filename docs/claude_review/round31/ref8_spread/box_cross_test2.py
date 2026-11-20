"""ILLUSTRATION (floating point, one random sample): prime-angle configurations of the VIII.3 fake
(angles concentrated), the VIII.4 spread fake (rational background + o(1/q) perturbation), Haar-random
angles, and the genuine angles of the same primes, at q = 383, P = q^1.99.
Tests:
 (i)  exact sector lattice count: for boxes [norm range of k consecutive cluster primes] x [angle interval],
      #model angles in box  <=  #Gaussian integers (one per associate class, both conjugates) in the sector.
      Genuine angles satisfy this trivially.  Reports the max ratio.
 (ii) max-scale cross gap at norm <= p: |sin(phi_j - arg w)| >= c / p_j for all w in Z[i], 0 < N w <= p_j,
      w not associate to the model's own factor direction.  Genuine angles satisfy it with c = 1
      (Im(pi_p conj w) is a nonzero integer).  Reports the number of primes failing at c = 1 and c = 0.1.
"""
import math, random, bisect, sys
import numpy as np
from sympy import isprime
sys.argv = [sys.argv[0]]
q = 383
rng = random.Random(11)
P = int(round(q ** 1.99))
pr = []
n = P
while len(pr) < q - 1:
    if n % 4 == 1 and isprime(n):
        pr.append(n)
    n += 1
print(f"P={P}, primes {pr[0]}..{pr[-1]}, sigma={math.log(pr[-1]/pr[0]):.4f}, sqrt(P)*sigma={math.sqrt(P)*math.log(pr[-1]/pr[0]):.1f}")
rng.shuffle(pr)  # random assignment of primes to labels 1..q-1
labels = list(range(1, q))
# ---- spread background (same algorithm as paley_spread.py, Astra's S choice) ----
sq = {x * x % q for x in range(1, q)}
chi = [0] + [1 if x in sq else -1 for x in range(1, q)]
K = np.array([[chi[(j - i) % q] for j in range(q)] for i in range(q)], dtype=np.int64)
A_ = [i for i in range(q) if K[i, 0] == 1]
S = A_[:-2]; Sset = set(S); Sc = [i for i in range(q) if i not in Sset]
from paley_spread import solve_mod_p
def near(nn):
    r = nn % (6 * q); r = min(r, 6 * q - r)
    return r ** 200 * q ** 199 <= (120 * q) ** 200
T = [j for j in range(q) if any(near(q + 4 * j + 4 * q * z) for z in range(3))]
tau1 = (K @ np.arange(q)) // q
u = np.zeros(q, dtype=np.int64); u[S] = tau1[S] % 3
x = solve_mod_p(K.T[np.ix_(T, Sc)], (-(K.T[np.ix_(T, S)] @ u[S])) % 3, 3)
u[Sc] = x[:, 0] % 3
z = (-pow(q, -1, 3) * (K.T @ u)) % 3
b = np.array([(q + 4 * j + 4 * q * int(z[j])) / (12 * q) for j in range(q)])
kappa, eps = 0.2, 0.01
c0 = rng.uniform(1 / 20, 1 / 10) / q ** 0.05
pert = np.array([(c0 + rng.uniform(-1, 1) * q ** (-kappa - eps) * rng.uniform(1, 2)) / (24 * q) for _ in range(q)])
theta_spread = (b + pert) % 1.0
H = q ** 0.08
theta_conc = np.array([(rng.uniform(1/(20*H), 1/(10*H)) * 0 + 0.07 + rng.uniform(-1, 1) * 2 * q ** (-0.243) / (100 * H * H)) for _ in range(q)])
theta_conc[:] = 0.07 + np.array([rng.uniform(-1, 1) for _ in range(q)]) * 2 * q ** (-0.243) / (100 * H * H)
theta_haar = np.array([rng.random() for _ in range(q)])
def gen_angle(p):
    for a in range(1, math.isqrt(p) + 1):
        bb = math.isqrt(p - a * a)
        if a * a + bb * bb == p:
            return math.atan2(bb, a) / (math.pi / 2)
theta_gen = np.zeros(q)
for j, p in zip(labels, pr):
    theta_gen[j] = gen_angle(p) if rng.random() < .5 else (1 - gen_angle(p)) % 1
models = {'genuine': theta_gen, 'concentrated(VIII.3)': theta_conc, 'spread(VIII.4)': theta_spread, 'Haar': theta_haar}
# ---- lattice points: one per associate class (a > 0, b >= 0), all norms in [P, max pr], angle in quarter turns
pmax = max(pr)
pts = []
for a in range(1, math.isqrt(pmax) + 1):
    lo = max(0, math.isqrt(max(0, P - a * a - 1)) - 1)
    for bb in range(lo, math.isqrt(max(0, pmax - a * a)) + 1):
        nn = a * a + bb * bb
        if P <= nn <= pmax:
            pts.append((nn, math.atan2(bb, a) / (math.pi / 2)))
pts.sort()
pnorm = np.array([t[0] for t in pts]); pang = np.array([t[1] for t in pts])
print("lattice classes in the annulus:", len(pts))
order = sorted(labels, key=lambda j: pr[j - 1])
for LMIN in (0, 5, 20):
  for name, th in models.items():
      worst = 0; worst_box = None
      for k in (1, 2, 4, 8, 16, 32, 64, 128, 381):
          for s0 in range(0, q - 1 - k + 1, max(1, k // 2)):
              sel = order[s0:s0 + k]
              n1, n2 = pr[sel[0] - 1], pr[sel[-1] - 1]
              mask = (pnorm >= n1) & (pnorm <= n2)
              la = np.sort(pang[mask])
              fa = np.sort(np.array([th[j] for j in sel]))
              # intervals spanned by consecutive model angles (all windows of r consecutive model angles)
              rs = sorted(set([1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, len(fa)]))
              for r in [x for x in rs if x <= len(fa)]:
                  for i0 in range(0, len(fa) - r + 1):
                      lo_, hi_ = fa[i0], fa[i0 + r - 1]
                      cnt_l = np.searchsorted(la, hi_ + 1e-9, 'right') - np.searchsorted(la, lo_ - 1e-9, 'left')
                      ratio = r / max(cnt_l, 0.5) if cnt_l >= LMIN else 0
                      if ratio > worst:
                          worst = ratio; worst_box = (k, r, cnt_l)
      print(f"(i) LMIN={LMIN} {name:22s} max (#model angles)/(#lattice classes) over boxes = {worst:.2f}  at (k primes, r angles, lattice)={worst_box}")
# ---- (ii) cross gap with w of norm <= p_j
allw = []
for a in range(-math.isqrt(pmax), math.isqrt(pmax) + 1):
    for bb in range(0, math.isqrt(pmax) + 1):
        nn = a * a + bb * bb
        if 0 < nn <= pmax and not (bb == 0 and a < 0):
            allw.append((math.atan2(bb, a), nn))
# reduce angles mod pi/2 (associates), keep norm
wa = np.array([(t[0] % (math.pi / 2)) for t in allw]); wn = np.array([t[1] for t in allw])
for name, th in models.items():
    fails = {1.0: 0, 0.1: 0}
    for j, p in zip(labels, pr):
        phi = th[j] * math.pi / 2
        m = wn <= p
        d = np.abs(np.sin(phi - wa[m]))
        # exclude the model's own direction class if genuine (distance ~0 to an associate of pi_p itself)
        own = (wn[m] == p) & (d < 1e-9)
        d = d[~own]
        dm = d.min()
        for c in fails:
            if dm < c / p:
                fails[c] += 1
    print(f"(ii) {name:22s} primes failing |sin(phi-arg w)| >= c/p for some N w <= p:  c=1: {fails[1.0]},  c=0.1: {fails[0.1]}  (of {q-1})")
