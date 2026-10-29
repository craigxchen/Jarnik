"""K8: random multi-cut slab boxes (round-11 convention, T = 1): criterion reach (first failure of SHIFT_e w.r.t.
L_pair, mod p) vs lattice-zonotope order zeta (exhaustive DFS, own code) vs rho.  Evidence only."""
import sys, random, itertools, math, time
sys.path.insert(0, '.')
from klib import build_box, Lpair, kfull
from klib3 import shift_dims
from klib import shiftset
def canon(n, S):
    Sc = tuple(x for x in range(n) if x not in S)
    return min(tuple(S), Sc, key=lambda s: (len(s), s))
def zeta(n, slabs, target, tlimit=600):
    d = n - 1
    def sval(v, S):
        vf = kfull(v); return sum(vf[x] for x in S)
    cands = []
    for v in itertools.product(range(-2, 3), repeat=d):
        if not any(v): continue
        g = 0
        for x in v: g = math.gcd(g, abs(x))
        if g != 1 or next(x for x in v if x) < 0: continue
        sw = [abs(sval(v, S)) for S, L in slabs]
        if any(w > 2 * L for w, (S, L) in zip(sw, slabs)): continue
        cands.append((v, sw))
    best = [0]; t0 = time.time(); timed = [False]
    def fits(gens):
        lo = [0] * d; hi = [0] * d
        for (v, sw), Lg in gens:
            for x in range(d):
                if v[x] > 0: hi[x] += Lg * v[x]
                else: lo[x] += Lg * v[x]
        rng = [range(-1 - lo[x], 1 - hi[x] + 1) for x in range(d)]
        if any(len(r) == 0 for r in rng): return False
        for k0 in itertools.product(*rng):
            ok = True
            for S, L in slabs:
                mn = sval(k0, S); mx = mn
                for (v, sw), Lg in gens:
                    s_ = sval(v, S) * Lg
                    if s_ > 0: mx += s_
                    else: mn += s_
                if mn < -L or mx > L: ok = False; break
            if ok: return True
        return False
    class Done(Exception): pass
    def dfs(start, gens, ccap, scap, tot):
        if time.time() - t0 > tlimit: timed[0] = True; raise Done
        if gens and not fits(gens): return
        if tot > best[0]:
            best[0] = tot
            if tot >= target: raise Done
        for i in range(start, len(cands)):
            v, sw = cands[i]
            for Lg in (1, 2):
                if any(Lg * abs(v[x]) > ccap[x] for x in range(d)): break
                if any(Lg * w > c for w, c in zip(sw, scap)): break
                dfs(i + 1, gens + [((v, sw), Lg)], [ccap[x] - Lg * abs(v[x]) for x in range(d)],
                    [c - Lg * w for w, c in zip(sw, scap)], tot + Lg)
    try:
        dfs(0, [], [2] * d, [2 * L for S, L in slabs], 0)
    except Done:
        pass
    return best[0], timed[0]
seed0 = int(sys.argv[1]); count = int(sys.argv[2]); n = int(sys.argv[3])
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
for sd in range(seed0, seed0 + count):
    rng = random.Random(sd)
    d = n - 1
    r = rng.randint(n - 1, 2 * n)
    pats = set()
    while len(pats) < r:
        S = tuple(x for x in range(n) if rng.random() < 0.5)
        if 0 < len(S) < n: pats.add(canon(n, S))
    L = rng.choice([1, 2])
    slabs = [(S, L) for S in sorted(pats)]
    A = build_box(n, slabs)
    H = []; reach = None
    for e in range(0, 2 * d + 2):
        Lr = [l for l in itertools.product((-1, 0, 1), repeat=d) if any(l) and len(shiftset(A, l)) > e]; hs, sm = shift_dims(A, Lr, e, 2147483629, d)
        H.append(hs)
        if hs == 0: break
        if reach is None and hs != sm: reach = e
    rho = len([h for h in H if h]) - 1
    if reach is None: reach = rho + 1
    z, timed = zeta(n, slabs, reach + 1)
    pr("seed=%d n=%d r=%d L=%d |A|=%d rho=%d H=%s rich-L criterion reach=%d zeta%s%d %s" % (
        sd, n, len(slabs), L, len(A), rho, H, reach, ">=" if (z > reach or timed) else "=", z,
        "(search timed out)" if timed else ("ZONOTOPE BEATS CRITERION" if z > reach else ("CRITERION BEATS ZONOTOPE" if z < reach else "equal"))))
open("k11_richrandom_%d_%d_%d_out.txt" % (seed0, count, n), "w").write("\n".join(out) + "\n")
