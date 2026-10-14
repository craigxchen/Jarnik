import itertools, time, numpy as np
from loclib import PR, image, pullback
from bigring import BigRing
from cyclib import CycRing, ideal_dims_cyc, hilbert_cyc
from jfast2 import jproj_gens_fast2
from syl16 import syl16_singer
import models
t0 = time.time(); m = 4
d, cols, Ls, colA = syl16_singer(1)
A = models.box(d, cols, Ls); ring = BigRing(d, m); cr = CycRing(ring, 15)
H = hilbert_cyc(cr, A, m, None); JA = [ring.N[i]-H[i] for i in range(m+1)]
print("Syl(16,1) |A|", len(A), "H", H, "J(A)", JA, flush=True)
coords = [[1 if y == x else 0 for y in range(d)] for x in range(d)]
Dset = cols + coords
def shift(i, t): return (i // 15) * 15 + (i % 15 + t) % 15   # columns are also shifted by the Singer cycle? check below
# verify column permutation under coordinate shift
colset = {tuple(c): i for i, c in enumerate(cols)}
perm = [colset[tuple(np.roll(c, 1))] for c in cols]
gens = {}
cum = lambda a: list(itertools.accumulate(a))
Av = np.array(A)
for b in range(1, 4):
    seen = set(); cnt = 0
    for S in itertools.combinations(range(30), b):
        if S in seen: continue
        orb = set(); cur = S
        for t in range(15):
            orb.add(tuple(sorted(cur)))
            cur = tuple((perm[i] if i < 15 else 15 + (i - 15 + 1) % 15) for i in cur)
        seen |= orb
        F = [Dset[i] for i in S]
        if np.linalg.matrix_rank(np.array(F)) < b: continue
        P = sorted(set(map(tuple, (Av @ np.array(F).T).tolist())))
        full = np.prod([len(set(p[z] for p in P)) for z in range(b)])
        if len(P) == full and b > 1: continue
        jp, rho, fails = jproj_gens_fast2(P, b, m); cnt += 1
        for i, forms in jp.items():
            for f in forms:
                vec = pullback(ring, F, f, i, PR, {})
                if vec.any(): gens.setdefault(i, []).append(vec)
    G = {k: np.array(v) for k, v in gens.items()}
    Jl = ideal_dims_cyc(cr, G, m)
    print(f"  mixed (cols+coords) b<={b} [{cnt} non-grid orbit reps]: {Jl} cum_eq={[x==y for x,y in zip(cum(Jl),cum(JA))]} [{time.time()-t0:.0f}s]", flush=True)
