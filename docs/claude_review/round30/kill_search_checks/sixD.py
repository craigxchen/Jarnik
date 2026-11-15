# kill_search: ALL 6-point rank-4 codes at D=10, c=0 (distance >= 5, balance-feasible), up to translation/perm symmetry.
import itertools, numpy as np, pickle, sys
sys.path.insert(0, '.')
import rankcodes as rc
from scipy.optimize import milp, LinearConstraint, Bounds
D = int(sys.argv[1]); L = D-1-D//2; h = -(-L//2); r = 4
zero = tuple([0]*r); basis = [tuple(1 if i == j else 0 for i in range(r)) for j in range(r)]
lams = [l for l in rc.Lam(r) if rc.feasible([zero]+basis+[l], r, D, L, h, True)]
print("feasible lambdas", len(lams))
def comps(n, k):
    # all k-tuples of nonneg ints summing to n
    for c in itertools.combinations(range(n+k-1), k-1):
        prev = -1; out = []
        for x in c: out.append(x-prev-1); prev = x
        out.append(n+k-1-prev-1); yield out
codes = set(); nbal = 0
def canon(P):  # P: 6 x D 0/1 matrix (points); symmetry: XOR-translation, coordinate perms, row perms
    best = None
    for t in range(6):
        Q = P ^ P[t]
        for perm in itertools.permutations(range(6)):
            key = tuple(sorted(tuple(col) for col in Q[list(perm)].T))
            if best is None or key < best: best = key
    return best
for lam in lams:
    A = rc.rowvals([zero]+basis+[lam], r)  # list of (u, values on 6 rows)
    vals = np.array([v for u, v in A])     # nA x 6
    nA = len(A)
    allm = np.array(list(comps(D, nA)), dtype=np.int16)
    ok = np.ones(len(allm), bool)
    pairs = list(itertools.combinations(range(6), 2))
    for i, j in pairs:
        dvec = (vals[:, i] != vals[:, j]).astype(np.int16)
        ok &= (allm @ dvec) >= L+1
    surv = allm[ok]
    for m in surv:
        # balance feasibility: split m_u = a_u + b_u (signs + / -)
        nv = 2*nA; Am = []; lo = []
        for i, j in pairs:
            a1 = np.zeros(nv); a2 = np.zeros(nv)
            for idx in range(nA):
                if vals[idx, i] < vals[idx, j]: a1[2*idx] = 1; a2[2*idx+1] = 1
                elif vals[idx, i] > vals[idx, j]: a2[2*idx] = 1; a1[2*idx+1] = 1
            Am += [a1, a2]; lo += [h, h]
        for idx in range(nA):
            e = np.zeros(nv); e[2*idx] = 1; e[2*idx+1] = 1; Am.append(e); lo.append(m[idx])
        hi = [np.inf]*(2*len(pairs)) + list(m)
        res = milp(c=np.zeros(nv), constraints=LinearConstraint(np.array(Am), lo, hi), integrality=np.ones(nv), bounds=Bounds(0, D))
        if res.status != 0: continue
        nbal += 1
        cols = []
        for idx in range(nA): cols += [vals[idx]]*int(m[idx])
        P = np.array(cols).T.astype(np.int8)  # 6 x D
        codes.add(canon(P))
    print("lambda", lam, "|A|", nA, "distance-ok", len(surv), "running balanced", nbal, "classes", len(codes), flush=True)
print("TOTAL 6-point rank-4 balanced-feasible codes up to symmetry:", len(codes))
pickle.dump(sorted(codes), open('six_D%d_codes.pkl' % D + '', 'wb'))
