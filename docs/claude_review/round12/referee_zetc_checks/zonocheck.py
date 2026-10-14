# Independent: rho_p (mod p rank profile, valid upper bound for rho_Q) and a MILP-found zonotope,
# verified point by point, for two-cut and three-cut blocks.
import itertools, sys, numpy as np, flint
from scipy.optimize import milp, LinearConstraint, Bounds
P = 2**61 - 1
def box(d, cuts):
    return [k for k in itertools.product((-1,0,1), repeat=d) if all(abs(sum(k[x] for x in S)) <= L for S, L in cuts)]
_mon = {}
def mons(d):
    if d not in _mon:
        m = [e for e in itertools.product((0,1,2), repeat=d)]; m.sort(key=sum); _mon[d] = m
    return _mon[d]
def rho_p(A, d):
    ms = mons(d)
    M = flint.nmod_mat(len(A), len(ms), [ (1 if all(b == 0 or a != 0 for a, b in zip(k, e)) and True else 0) and
        int(np.prod([a**b for a, b in zip(k, e)])) % P for k in A for e in ms], P)
    R, r = M.rref()
    assert r == len(A)
    last = 0
    for i in range(r):
        for j in range(len(ms)):
            if R[i, j] != 0:
                last = max(last, sum(ms[j])); break
    return last
def zono(d, cuts, A):
    cand = []
    for v in itertools.product((-1,0,1), repeat=d):
        if any(v) and next(x for x in v if x) > 0: cand.append(v)
    funcs = [tuple(1 if y == x else 0 for y in range(d)) for x in range(d)]
    bnds = [1]*d
    for S, L in cuts:
        funcs.append(tuple(1 if y in S else 0 for y in range(d))); bnds.append(L)
    nv = len(cand); nvar = nv + d
    rows = []; ub = []
    for f, b in zip(funcs, bnds):
        fv = [sum(a*c for a, c in zip(f, v)) for v in cand]
        # upper: f(k0) + sum L_v max(0,fv) <= b
        rows.append([max(0, x) for x in fv] + list(f)); ub.append(b)
        # lower: -f(k0) + sum L_v max(0,-fv) <= b
        rows.append([max(0, -x) for x in fv] + [-a for a in f]); ub.append(b)
    c = np.array([-1.0]*nv + [0.0]*d)
    res = milp(c, constraints=LinearConstraint(np.array(rows, float), -np.inf, np.array(ub, float)),
               integrality=np.ones(nvar), bounds=Bounds([0]*nv + [-1]*d, [4]*nv + [1]*d))
    x = np.round(res.x).astype(int)
    edges = [(cand[i], x[i]) for i in range(nv) if x[i] > 0]; k0 = tuple(x[nv:])
    Aset = set(A)
    # verify all lattice points
    for ts in itertools.product(*[range(Lv + 1) for _, Lv in edges]):
        pt = tuple(k0[j] + sum(t * v[j] for t, (v, _) in zip(ts, edges)) for j in range(d))
        assert pt in Aset, "zonotope point outside A"
    return sum(Lv for _, Lv in edges), edges, k0
