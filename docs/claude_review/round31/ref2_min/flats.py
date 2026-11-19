"""Independent BFS census of relation flats (referee code, ref2_min).
Rows: e_i+e_l-e_j-e_k (i<j<k<l), optionally midpoint rows e_i+e_k-2e_j (i<j<k).
Every flat of rank k+1 is cl(F + r) for an admissible rank-k flat F (admissibility is
inherited by subflats), so BFS over closures with dedupe is exhaustive.
Arithmetic mod P=10007: all ranks are ranks of integer vectors with Euclidean norm <= sqrt6,
sum zero, so nonzero minors have |det| <= 6^(7/2) < 300 < P; ranks mod P equal ranks over Q.
Counts reported: n70 (generating rows in flat), n210 (all 4-index pairings), nmid (56 midpoint rows
in all orientations: 2e_j - e_i - e_k for every j and pair {i,k} not containing j: 8*21=168).
"""
import sys, itertools, numpy as np, pickle
P = 10007
mid = len(sys.argv) > 1 and sys.argv[1] == 'mid'
def vec(d):
    v = np.zeros(8, dtype=np.int64)
    for k, c in d.items(): v[k] += c
    return v % P
gen = [vec({i:1,l:1,j:-1,k:-1}) for i,j,k,l in itertools.combinations(range(8),4)]
if mid:
    gen += [vec({i:1,k:1,j:-2}) for i,j,k in itertools.combinations(range(8),3)]
GEN = np.array(gen)
ALL4 = []
for q in itertools.combinations(range(8),4):
    a,b,c,d = q
    for (p1,p2) in [((a,b),(c,d)),((a,c),(b,d)),((a,d),(b,c))]:
        ALL4.append(vec({p1[0]:1,p1[1]:1,p2[0]:-1,p2[1]:-1}))
ALL4 = np.array(ALL4)
MIDS = np.array([vec({i:1,k:1,j:-2}) for j in range(8) for i,k in itertools.combinations([x for x in range(8) if x!=j],2)])
INV = [0]+[pow(x,P-2,P) for x in range(1,P)]

def reduce_rows(M, basis):  # basis: list of (pivot, row) RREF-like (each row has 1 at pivot, 0 at other pivots)
    M = M.copy()
    for piv, row in basis:
        M = (M - np.outer(M[:,piv], row)) % P
    return M

def add(basis, v):
    v = reduce_rows(v[None,:], basis)[0]
    nz = np.nonzero(v)[0]
    if len(nz) == 0: return None
    piv = nz[0]; v = v*INV[v[piv]] % P
    nb = []
    for p2, r in basis:
        r = (r - r[piv]*v) % P
        nb.append((p2, r))
    nb.append((piv, v))
    return nb

QUADS = list(itertools.combinations(range(8),4))
def admissible(basis):
    n = reduce_rows(np.eye(8, dtype=np.int64), basis)  # normal forms of e_i
    for i,j in itertools.combinations(range(8),2):
        if not ((n[i]-n[j]) % P).any(): return False
    for i,j,k,l in QUADS:
        D = np.array([n[j]-n[i], n[k]-n[i], n[l]-n[i]]) % P
        # rank <=1 iff all 2x2 minors vanish
        ok = True
        for a,b in ((0,1),(0,2),(1,2)):
            m = (np.outer(D[a],D[b]) - np.outer(D[b],D[a])) % P
            if m.any(): ok = False; break
        if ok: return False
    return True

def counts(basis):
    return (int((~reduce_rows(GEN,basis).any(1)).sum()), int((~reduce_rows(ALL4,basis).any(1)).sum()),
            int((~reduce_rows(MIDS,basis).any(1)).sum()))

def key(basis):
    z = ~reduce_rows(GEN, basis).any(1)
    return sum(1<<int(i) for i in np.nonzero(z)[0])

level = {0: []}  # key -> basis
census = {}
for rank in range(1, 6):
    new = {}
    for k0, basis in level.items():
        for r in range(len(GEN)):
            if (k0 >> r) & 1: continue
            nb = add(basis, GEN[r])
            kk = key(nb)
            if kk in new: continue
            new[kk] = nb if admissible(nb) else None
    level = {k: b for k, b in new.items() if b is not None}
    hist = {}
    recs = []
    for k, b in level.items():
        c = counts(b)
        hist[c] = hist.get(c, 0) + 1
        recs.append((k, c, [(int(p), [int(x) for x in r]) for p, r in b]))
    census[rank] = recs
    print('rank', rank, 'flats total', len(new), 'admissible', len(level), flush=True)
    print('  max n70', max(c[0] for c in hist), 'max n210', max(c[1] for c in hist), 'max n4+mid', max(c[1]+c[2] for c in hist))
    print('  hist (n70,n210,nmid)->#:', sorted(hist.items(), key=lambda t: -t[0][1]-t[0][2])[:12], flush=True)
pickle.dump(census, open('census_%s.pkl' % ('mid' if mid else 'tern'), 'wb'))
