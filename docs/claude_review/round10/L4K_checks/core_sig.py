# Exhaustive (level 1) enumeration of ALL good systems on a set of small forced moduli ("core"), with budgets
# B(v,s) for all ||v||_1 <= K, and classification: at each q, is the label F1-equivalent to genuine
# (c = u L mod n, u unit)?  Usage: python3 core_sig.py K r tag q1,q2,...
import sys, math, time
import numpy as np
from lib import *
from forced_sig import sig_table, budgets
from forced_coh import load_set

K, r, tag = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
core = [int(x) for x in sys.argv[4].split(',')]
P, F = load_set(K, r, tag)
t0 = time.time()
tables = [sig_table(P, K, q) for q in core]
V = tables[0]['V']; B = budgets(P, V); nv = len(V)
rest = [q for q in F if q not in core]
if len(sys.argv) > 5 and sys.argv[5] == 'rest':
    # subtract the genuine all-level demand of the forced moduli outside the core (they stay genuine)
    def vq(x, q):
        x = abs(x); k = 0
        while x and x % q == 0: x //= q; k += 1
        return k
    pis = [gauss_prime(p) for p in P]
    for i, v in enumerate(V):
        z = Zv(pis, v)
        for s in range(4):
            B[i, s] -= sum(vq(Ds(z, s), q) * math.log(q) for q in rest)
    print(f"  budgets reduced by the genuine demand of the other forced moduli {rest}")
if len(sys.argv) > 6:   # replace the genuine demand of the saved top moduli by the saved non-genuine system's load
    import pickle
    pk = pickle.load(open(sys.argv[6], 'rb')); sv = pk['saved'][int(sys.argv[7])]
    top = sorted(sv['system'])
    def vq2(x, q):
        x = abs(x); k = 0
        while x and x % q == 0: x //= q; k += 1
        return k
    pis2 = [gauss_prime(p) for p in P]
    for i, v in enumerate(V):
        z = Zv(pis2, v)
        for s in range(4):
            B[i, s] += sum(vq2(Ds(z, s), q) * math.log(q) for q in top)   # undo genuine top demand
    B -= sv['top_load']
    print(f"  top moduli {top} replaced by saved non-genuine system; min remaining budget {B.min():.3f}")
# per q: which signature classes contain only genuine-class vectors / some
genclass_keys = []
for t in tables:
    n = t['n']; Lc = t['L'] % n
    units = [u for u in range(1, n + 1) if math.gcd(u, n) == 1] if n > 1 else [1]
    G = set(tuple(int(x) for x in (u * Lc) % n) for u in units)
    isg = np.array([tuple(int(x) for x in row) in G for row in (t['ells'] % n)])
    # per key: number of genuine-class vectors and total
    ng = np.bincount(t['inv'], weights=isg, minlength=len(t['keys']))
    genclass_keys.append(ng)
costs = []
for t in tables:
    lq = math.log(t['q']); C = []
    for key in t['keys']:
        c = np.zeros((nv, 4))
        for i, s in enumerate(key):
            if s >= 0: c[i, s] = lq
        C.append(c)
    costs.append(C)
order = list(range(len(tables)))
tot = [0]; totgen = [0]; ntup = [0]; per_q_nongen = [0] * len(tables)
def dfs(k, load, mult, multgen, choice):
    if k == len(order):
        tot[0] += mult; totgen[0] += multgen; ntup[0] += 1
        for i, kk in choice.items():
            t = tables[i]; cnt = t['counts'][kk]; g = genclass_keys[i][kk]
            # systems whose label at q_i is non-genuine-class: mult * (cnt - g)/cnt
            per_q_nongen[i] += mult * (cnt - g) / cnt
        return
    i = order[k]; t = tables[i]
    for kk, c in enumerate(costs[i]):
        nl = load + c
        if (nl <= B + 1e-9).all():
            choice[i] = kk
            dfs(k + 1, nl, mult * int(t['counts'][kk]), multgen * int(genclass_keys[i][kk]), choice)
            del choice[i]
dfs(0, np.zeros((nv, 4)), 1, 1, {})
print(f"K={K} r={r} {tag} P={P} core={core}: good level-1 systems on core {tot[0]} ({ntup[0]} signature tuples); "
      f"genuine class at every core q: {totgen[0]}; [{time.time()-t0:.0f}s]")
for i, t in enumerate(tables):
    print(f"   q={t['q']}: good systems with non-genuine-class label at q: {per_q_nongen[i]:.0f}")
