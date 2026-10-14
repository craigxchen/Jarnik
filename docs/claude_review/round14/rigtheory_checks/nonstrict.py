"""Which genuine-class labels survive one-modulus deformation at q0 (exact): ell' = u*L + n*w mod m."""
import sys, math, itertools
from rlib14 import *
from forced import forced_set
from deform import setup
K = int(sys.argv[1]); ps = [int(x) for x in sys.argv[2].split(",")]; q0 = int(sys.argv[3])
F = forced_set(ps, K, 60000)
V, m, L, alw, gcoll, gsig, last = setup(ps, K, F, q0)
n = m // 4; r = len(ps); nv = nu(ps, q0)
surv = {}
for u in range(1, m):
    if math.gcd(u, m) != 1: continue
    for w in itertools.product(range(4), repeat=r):
        ell = [(u * L[j] + n * w[j]) % m for j in range(r)]
        if any(ell[j] % 2 != nv[j] for j in range(r)): continue
        ok = True
        for i, v in enumerate(V):
            t = sum(a * b for a, b in zip(v, ell)) % m
            if t % n == 0 and not alw[i, t // n]: ok = False; break
        if ok:
            key = tuple(ell)
            surv.setdefault(key, set()).add((u % 4, w))
print(f"K={K} q0={q0} n={n}: surviving class labels={len(surv)}")
types = {}
for ell, reps in surv.items():
    # does some representation have w == 0 ?
    t = ("u=1mod4,w=0" if any(a == 1 and all(x == 0 for x in w) for a, w in reps) else
         "u=3mod4,w=0" if any(a == 3 and all(x == 0 for x in w) for a, w in reps) else "needs w!=0")
    types[t] = types.get(t, 0) + 1
print(types)
