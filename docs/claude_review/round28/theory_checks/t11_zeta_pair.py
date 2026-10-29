"""T11: zeta_pair(A) = max order of a lattice zonotope k0 + sum t_i v_i (t_i in {0..m_i}, m_i in {1,2}) inside A with
ALL generators in L_pair (the criterion's own hypothesis class).  Exhaustive DFS (fitting is monotone under adding
generators), adapted from the round-27 killer's k6_zono.py with candidate set restricted to L_pair.
Usage: t11_zeta_pair.py model L   (model syl8 | sdc7)"""
import sys, itertools, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
model, L = sys.argv[1], int(sys.argv[2])
n, pats = syl(3) if model == 'syl8' else sdc(7)
slabs = [(S, L) for S in pats]
d = n - 1
A = set(build_box(n, slabs))
def sval(v, S):
    vf = kfull(v); return sum(vf[x] for x in S)
cands = [(v, [abs(sval(v, S)) for S, _ in slabs]) for v in Lpair_half(d)]
best = [0, None]; t0 = time.time(); nodes = [0]
def fits(gens):
    lo = [0] * d; hi = [0] * d
    for (v, sw), m in gens:
        for x in range(d):
            if v[x] > 0: hi[x] += m * v[x]
            else: lo[x] += m * v[x]
    rng = [range(-1 - lo[x], 1 - hi[x] + 1) for x in range(d)]
    if any(len(r) == 0 for r in rng): return None
    for k0 in itertools.product(*rng):
        ok = True
        for S, LL in slabs:
            mn = sval(k0, S); mx = mn
            for (v, sw), m in gens:
                s_ = sval(v, S) * m
                if s_ > 0: mx += s_
                else: mn += s_
            if mn < -LL or mx > LL: ok = False; break
        if ok: return k0
    return None
def dfs(start, gens, ccap, scap, tot):
    nodes[0] += 1
    if gens:
        k0 = fits(gens)
        if k0 is None: return
        if tot > best[0]:
            best[0] = tot; best[1] = (k0, [(g[0], m) for g, m in gens])
    for i in range(start, len(cands)):
        v, sw = cands[i]
        for m in (1, 2):
            if any(m * abs(v[x]) > ccap[x] for x in range(d)): break
            if any(m * w > c for w, c in zip(sw, scap)): break
            dfs(i + 1, gens + [((v, sw), m)], [ccap[x] - m * abs(v[x]) for x in range(d)],
                [c - m * w for w, c in zip(sw, scap)], tot + m)
dfs(0, [], [2] * d, [2 * LL for S, LL in slabs], 0)
print("%s L=%d |A|=%d: zeta_pair = %d (exhaustive over L_pair generators, multiplicity <= 2; nodes %d; %.0fs)" % (
      model, L, len(A), best[0], nodes[0], time.time() - t0))
print("  witness k0 =", best[1][0], " generators:", best[1][1])
