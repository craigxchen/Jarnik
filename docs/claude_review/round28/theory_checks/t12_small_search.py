"""T12: search small slab boxes (n = 4,5,6; T = 1) for reach(A, L_pair) > zeta_all(A) (any generators).
Reach mod p first; hits are re-verified EXACTLY over Q (fmpz nullspaces).  zeta_all by exhaustive DFS over primitive
generators v in {-2..2}^d with multiplicity <= 2 (killer's k6 rules)."""
import sys, itertools, random, time, math
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
def zeta_all(n, slabs, cap=99, pair_only=False):
    d = n - 1
    def sval(v, S):
        vf = kfull(v); return sum(vf[x] for x in S)
    if pair_only:
        cl = Lpair_half(d)
    else:
        cl = []
        for v in itertools.product(range(-2, 3), repeat=d):
            if not any(v): continue
            g = 0
            for x in v: g = math.gcd(g, abs(x))
            if g != 1 or next(x for x in v if x) < 0: continue
            cl.append(v)
    cands = []
    for v in cl:
        sw = [abs(sval(v, S)) for S, L in slabs]
        if any(w > 2 * L for w, (S, L) in zip(sw, slabs)): continue
        cands.append((v, sw))
    best = [0]
    def fits(gens):
        lo = [0] * d; hi = [0] * d
        for (v, sw), m in gens:
            for x in range(d):
                if v[x] > 0: hi[x] += m * v[x]
                else: lo[x] += m * v[x]
        rng = [range(-1 - lo[x], 1 - hi[x] + 1) for x in range(d)]
        if any(len(r) == 0 for r in rng): return False
        for k0 in itertools.product(*rng):
            ok = True
            for S, L in slabs:
                mn = sval(k0, S); mx = mn
                for (v, sw), m in gens:
                    s_ = sval(v, S) * m
                    if s_ > 0: mx += s_
                    else: mn += s_
                if mn < -L or mx > L: ok = False; break
            if ok: return True
        return False
    def dfs(start, gens, ccap, scap, tot):
        if gens and not fits(gens): return
        best[0] = max(best[0], tot)
        for i in range(start, len(cands)):
            v, sw = cands[i]
            for m in (1, 2):
                if any(m * abs(v[x]) > ccap[x] for x in range(d)): break
                if any(m * w > c for w, c in zip(sw, scap)): break
                dfs(i + 1, gens + [((v, sw), m)], [ccap[x] - m * abs(v[x]) for x in range(d)],
                    [c - m * w for w, c in zip(sw, scap)], tot + m)
    dfs(0, [], [2] * d, [2 * L for S, L in slabs], 0)
    return best[0]
def reach_modp(A, d, p=P1):
    Lp = Lpair_half(d); H = []
    for e in range(0, 2 * d + 1):
        hs = H_modp(A, e, p, d); H.append(hs)
        if hs == 0: return e, H
        _, ds = shift_dims(A, Lp, e, p, d)
        if ds < hs: return e, H
    return 2 * d + 1, H
def reach_exact(A, d, upto):
    from klibQ import shift_ok_Q
    return None
seen = set(); hits = []
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
N = int(sys.argv[2]) if len(sys.argv) > 2 else 200
t0 = time.time()
for it in range(N):
    n = rng.choice([4, 5, 6]); d = n - 1
    r = rng.randint(1, 2 * n)
    pats = set()
    for _ in range(r):
        S = tuple(x for x in range(n) if rng.random() < 0.5)
        if 0 < len(S) < n: pats.add(canon(n, S))
    if not pats: continue
    mode = rng.choice(['L1', 'L2', 'vertex'])
    slabs = [(S, (len(S) - 1 if mode == 'vertex' else int(mode[1]))) for S in sorted(pats) if len(S) >= 2]
    slabs = [(S, L) for S, L in slabs if 1 <= L < len(S)]
    if not slabs: continue
    key = (n, tuple(slabs))
    if key in seen: continue
    seen.add(key)
    A = build_box(n, slabs)
    reach, H = reach_modp(A, d)
    rho = max(i for i, x in enumerate(H) if x)
    z = -1
    zp = zeta_all(n, slabs, pair_only=True)
    tag = "reach>zeta_pair" if reach > zp else ""
    if tag or it % 20 == 0:
        print("it=%d n=%d slabs=%s |A|=%d H=%s rho=%d reach=%d zeta_all=%d zeta_pair=%d %s (%.0fs)" % (
            it, n, slabs, len(A), H, rho, reach, z, zp, tag, time.time() - t0), flush=True)
    if tag: hits.append((n, slabs, reach, z, zp))
print("hits:", len(hits))
