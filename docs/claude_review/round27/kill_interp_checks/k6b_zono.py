"""K6: lattice-zonotope order zeta(A) (round 11 Prop 2.1 sense: k0 + sum t_i v_i, 0<=t_i<=L_i, all in A)
for slab boxes in the round-11 convention.  Exhaustive DFS over generators v in {-2..2}^d (k'-coords),
primitive, with per-coordinate and per-slab width budgets; existence of k0 checked by brute force.
Usage: k6_zono.py model target   (model in syl8, sdc7, cut:d,s,L)."""
import sys, itertools, math, time
sys.path.insert(0, '.')
from klib import build_box, kfull
def canon(n, S):
    Sc = tuple(x for x in range(n) if x not in S)
    return min(tuple(S), Sc, key=lambda s: (len(s), s))
def syl8():
    H = [[1]]
    for _ in range(3):
        H = [r + r for r in H] + [r + [-x for x in r] for r in H]
    return tuple(sorted(set(canon(8, tuple(x for x in range(8) if r[x] == 1)) for r in H[1:])))
import numpy as np
def legendre(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)
def sdc7():
    q = 7; idx = ['inf'] + list(range(q)); M = q + 1
    K = np.zeros((M, M), dtype=int)
    for a, x in enumerate(idx):
        for b, y in enumerate(idx):
            if x == 'inf': K[a, b] = 1
            elif y == 'inf': K[a, b] = -1
            else: K[a, b] = (1 if x == y else 0) + legendre(y - x, q)
    if not ((K + K.T) == 2 * np.eye(M, dtype=int)).all(): K = K.T
    Dg = np.diag(K[:, 0]); K = Dg @ K @ Dg
    pats = set()
    for a in range(1, M):
        S = [x for x in range(M) if K[x, a] == -1]
        pats.add(canon(M, tuple(S))); pats.add(canon(M, tuple(sorted(set(S) ^ {0, a}))))
    return tuple(sorted(pats))
model = sys.argv[1]; target = int(sys.argv[2])
if model == 'syl8': n, slabs = 8, [(S, 2) for S in syl8()]
elif model == 'syl8L3': n, slabs = 8, [(S, 3) for S in syl8()]
elif model == 'sdc7': n, slabs = 8, [(S, 2) for S in sdc7()]
else:
    dd, s, L = (int(x) for x in model.split(':')[1].split(',')); n = dd + 1; slabs = [(tuple(range(s)), L)]
d = n - 1
A = set(build_box(n, slabs))
def sval(v, S):
    vf = kfull(v); return sum(vf[x] for x in S)
cands = []
for v in itertools.product(range(-2, 3), repeat=d):
    if not any(v): continue
    g = 0
    for x in v: g = math.gcd(g, abs(x))
    if g != 1: continue
    # canonical sign: first nonzero positive
    fnz = next(x for x in v if x)
    if fnz < 0: continue
    sw = [abs(sval(v, S)) for S, L in slabs]
    if any(w > 2 * L for w, (S, L) in zip(sw, slabs)): continue
    cands.append((v, sw))
print("model", model, "|A|", len(A), "candidate generators", len(cands), flush=True)
coordcap = [2] * d; slabcap = [2 * L for S, L in slabs]
best = [0, None]
def fits(gens):
    # does some k0 exist with all zonotope lattice points in A?  range per coordinate and slab
    lo = [0] * d; hi = [0] * d
    for (v, sw), Lg in gens:
        for x in range(d):
            if v[x] > 0: hi[x] += Lg * v[x]
            else: lo[x] += Lg * v[x]
    rng = [range(-1 - lo[x], 1 - hi[x] + 1) for x in range(d)]
    if any(len(r) == 0 for r in rng): return None
    for k0 in itertools.product(*rng):
        ok = True
        for S, L in slabs:
            mn = sval(k0, S); mx = mn
            for (v, sw), Lg in gens:
                s_ = sval(v, S) * Lg
                if s_ > 0: mx += s_
                else: mn += s_
            if mn < -L or mx > L: ok = False; break
        if ok: return k0
    return None
t0 = time.time(); nodes = [0]
def dfs(start, gens, ccap, scap, tot):
    nodes[0] += 1
    if tot > best[0]:
        k0 = fits(gens)
        if k0 is not None:
            best[0] = tot; best[1] = (k0, [(g[0], Lg) for g, Lg in gens])
            print("  zonotope of order", tot, ":", best[1], "(%.0fs)" % (time.time() - t0), flush=True)
            if tot >= target: raise StopIteration
        else:
            return   # fitting is monotone: if gens does not fit, no superset fits
    for i in range(start, len(cands)):
        v, sw = cands[i]
        for Lg in (1, 2):
            if any(Lg * abs(v[x]) > ccap[x] for x in range(d)): break
            if any(Lg * w > c for w, c in zip(sw, scap)): break
            dfs(i + 1, gens + [((v, sw), Lg)], [ccap[x] - Lg * abs(v[x]) for x in range(d)],
                [c - Lg * w for w, c in zip(sw, scap)], tot + Lg)
try:
    dfs(0, [], coordcap, slabcap, 0)
except StopIteration:
    pass
print("model %s: max zonotope order found %d (target %d); exhaustive=%s; nodes %d; %.0fs" % (
    model, best[0], target, best[0] < target, nodes[0], time.time() - t0), flush=True)
