import numpy as np, itertools
from rlib import chi_table
out = []
def log(s): print(s, flush=True); out.append(s)
def circ(q):
    ch = chi_table(q)
    return np.array([[ch[(a-t) % q] for t in range(q)] for a in range(q)], dtype=np.int64)  # C[a,t]=chi(a-t)
def j0(f, q):
    T = np.nonzero(f)[0]
    for j in range(q):
        if sum(int(f[t])*pow(int(t), j, q) for t in T) % q: return j
    return None
# exhaustive small q
for q, vals in [(7, range(-3, 4)), (11, range(-1, 2)), (19, range(-1, 2))]:
    C = circ(q); vals = list(vals); worst = -10**9; tight = 0; cnt = 0; viol = 0
    if q == 19:
        # ternary with support <= 10, sampled exhaustively over supports of size 4 and random larger
        rng = np.random.default_rng(1); F = []
        for _ in range(200000):
            k = rng.integers(1, (q+1)//2 + 1); f = np.zeros(q, dtype=np.int64)
            f[rng.choice(q, k, replace=False)] = rng.choice([-1, 1], k); F.append(f)
        F = np.array(F)
    else:
        F = np.array(list(itertools.product(vals, repeat=q)), dtype=np.int64)
        F = F[((F != 0).sum(1) >= 1) & ((F != 0).sum(1) <= (q+1)//2)]
    Z = ((F @ C.T) == 0).sum(1)
    for f, z in zip(F, Z):
        if not (f % q != 0).any(): continue
        jj = j0(f, q); bound = (q-1)//2 - jj; cnt += 1
        if z > bound: viol += 1
        if z == bound: tight += 1
        worst = max(worst, z - bound)
    log(f"q={q} vals={vals[0]}..{vals[-1]} : {cnt} f tested, violations of #zeros <= (q-1)/2 - j0: {viol}, max(z-bound)={worst}, tight cases={tight}")
# large values incl. multiples of q
rng = np.random.default_rng(2)
for q in [23, 31, 43, 59, 83, 103]:
    C = circ(q); viol = 0; worst = -10**9
    for _ in range(3000):
        k = rng.integers(1, (q+1)//2 + 1); f = np.zeros(q, dtype=np.int64)
        f[rng.choice(q, k, replace=False)] = rng.integers(-3*q, 3*q+1, k)
        if rng.random() < 0.3: f[f != 0][:1] = q
        if not (f % q != 0).any(): continue
        z = int(((C @ f) == 0).sum()); b = (q-1)//2 - j0(f, q); worst = max(worst, z-b); viol += z > b
    log(f"q={q}: 3000 random f (values up to 3q), violations {viol}, max(z-bound) {worst}")
# does the support condition matter? supports > (q+1)/2 with F=0: max zeros
for q in [11, 19, 23]:
    C = circ(q); best = 0; rng = np.random.default_rng(3)
    for _ in range(100000):
        k = rng.integers((q+3)//2, q+1); f = np.zeros(q, dtype=np.int64); f[rng.choice(q, k, replace=False)] = rng.choice([-1, 1], k)
        if f.sum() != 0: continue
        best = max(best, int(((C @ f) == 0).sum()))
    log(f"q={q}: F=0 ternary, support > (q+1)/2: max #zeros found {best} vs (q-3)/2 = {(q-3)//2}")
open('check3_out.txt', 'w').write('\n'.join(out) + '\n')
