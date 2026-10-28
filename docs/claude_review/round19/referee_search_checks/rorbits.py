"""Independent S_7-orbit representatives of D-multisets of edges of K_7 (numpy min-code canonical form)
and a Burnside cross-check of the orbit count.  Usage: python3 rorbits.py D"""
import sys, json, itertools, math
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/referee_search_checks')
from rcore import EDGES, EIDX
D = int(sys.argv[1])
perms = list(itertools.permutations(range(7)))
# sigma acts on edges
EP = np.array([[EIDX[tuple(sorted((s[i], s[j])))] for (i, j) in EDGES] for s in perms], dtype=np.int64)  # 5040x21
POW = np.array([7 ** k for k in range(21)], dtype=np.int64)
W = POW[EP].T  # 21 x 5040 : W[e, s] = 7^{sigma(e)}
mults = list(itertools.combinations_with_replacement(range(21), D))
V = np.zeros((len(mults), 21), dtype=np.int64)
for r, m in enumerate(mults):
    for e in m:
        V[r, e] += 1
canon = np.empty(len(mults), dtype=np.int64)
CH = 4000
for s in range(0, len(mults), CH):
    canon[s:s + CH] = (V[s:s + CH] @ W).min(axis=1)
u, first = np.unique(canon, return_index=True)
reps = [[list(EDGES[e]) for e in mults[i]] for i in first]
# Burnside
def fixed(cyc):
    # coefficient of t^D in prod 1/(1-t^l)
    poly = [1] + [0] * D
    for l in cyc:
        for d in range(l, D + 1):
            poly[d] += poly[d - l]
    return poly[D]
tot = 0
for row in EP:
    seen = [False] * 21; cyc = []
    for e in range(21):
        if not seen[e]:
            l = 0; f = e
            while not seen[f]:
                seen[f] = True; f = row[f]; l += 1
            cyc.append(l)
    tot += fixed(cyc)
assert tot % 5040 == 0
npairs = 0
for i in first:
    s = int(np.count_nonzero(V[i]))
    npairs += math.comb(21 - s + D - 1, D)
print(json.dumps(dict(D=D, multisets=len(mults), orbits=len(reps), burnside=tot // 5040, pairs=npairs)))
json.dump(reps, open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/referee_search_checks/rreps_D%d.json' % D, 'w'))
