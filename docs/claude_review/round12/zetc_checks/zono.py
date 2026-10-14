"""Zonotope lower bounds for zeta(A) on T = 1 slab-boxes, with explicit verification.

A = {k in {-1,0,1}^d : |k(S_p)| <= L_p}.  A lattice zonotope Z = k0 + sum_i [0, l_i] v_i (integers l_i >= 1)
lies in A iff every point does; we search by MILP over multiplicities l_v of primitive directions v in A - A with
the necessary width constraints  sum_v l_v |<v, f>| <= 2 c_f  (f = coordinates (c=1) and columns (c = L_p)),
then look for a base point k0 with all points of Z in A (exhaustive over k0 in {-1,0,1}^d, points enumerated
explicitly).  If placement fails, the MILP solution is cut off and re-solved (no-good cut).  Every reported value
is backed by an explicitly verified zonotope (all its lattice points checked), so it is an exact lower bound.
"""
import itertools, math
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds


def prim_dirs(pts):
    P = np.array(pts, dtype=np.int64)
    S = set()
    for i in range(len(P)):
        D = P[i] - P[i + 1:]
        for v in D:
            g = 0
            for x in v:
                g = math.gcd(g, int(abs(x)))
            if g == 0:
                continue
            v = tuple(int(x) // g for x in v)
            for x in v:
                if x != 0:
                    if x < 0:
                        v = tuple(-y for y in v)
                    break
            S.add(v)
    return sorted(S)


def zonotope_points(k0, edges):
    """edges: list of (v, l). Returns set of all lattice points k0 + sum t_i v_i, 0<=t_i<=l_i."""
    pts = {tuple(k0)}
    for v, l in edges:
        new = set()
        for p in pts:
            for t in range(l + 1):
                new.add(tuple(p[j] + t * v[j] for j in range(len(p))))
        pts = new
    return pts


def place(edges, d, Aset):
    """find k0 with zonotope inside A; return k0 or None."""
    # necessary: per coordinate, spread sum l|v_x| <= 2; choose k0_x so that range within [-1,1]
    lo = [0] * d; hi = [0] * d
    for v, l in edges:
        for j in range(d):
            if v[j] > 0: hi[j] += l * v[j]
            else: lo[j] += l * v[j]
    ranges = []
    for j in range(d):
        # k0_j + lo_j >= -1, k0_j + hi_j <= 1
        a = -1 - lo[j]; b = 1 - hi[j]
        if a > b:
            return None
        ranges.append(range(a, b + 1))
    for k0 in itertools.product(*ranges):
        Z = zonotope_points(k0, edges)
        if all(z in Aset for z in Z):
            return k0
    return None


def zeta_lower(pts, cols, Ls, d, dirs=None, avoid=None, max_iter=60, time_limit=60):
    """Max verified zonotope degree. avoid: optional predicate v -> bool; directions with avoid(v) True are
    excluded (used for nu lower bounds at a point u: exclude v with <v,u> = 0)."""
    Aset = set(pts)
    if dirs is None:
        dirs = prim_dirs(pts)
    if avoid is not None:
        dirs = [v for v in dirs if not avoid(v)]
    if not dirs:
        return 0, []
    funcs = [tuple(1 if x == y else 0 for x in range(d)) for y in range(d)]
    caps = [2] * d
    for S, L in zip(cols, Ls):
        funcs.append(tuple(1 if x in S else 0 for x in range(d)))
        caps.append(2 * L)
    F = np.array(funcs, dtype=float)
    V = np.array(dirs, dtype=float)
    M = np.abs(F @ V.T)
    cons = [LinearConstraint(M, -np.inf, np.array(caps, dtype=float))]
    nogoods = []
    best = (0, [])
    for it in range(max_iter):
        allc = list(cons)
        for ng in nogoods:
            allc.append(ng)
        res = milp(c=-np.ones(len(dirs)), constraints=allc, integrality=np.ones(len(dirs)),
                   bounds=Bounds(0, 2), options={'time_limit': time_limit})
        if res.x is None:
            break
        x = np.round(res.x).astype(int)
        val = int(x.sum())
        if val <= best[0]:
            break
        edges = [(dirs[i], int(x[i])) for i in range(len(dirs)) if x[i] > 0]
        k0 = place(edges, d, Aset)
        if k0 is not None:
            best = (val, edges)
            break
        # no-good cut: forbid this exact multiset (sum over its support of x_i <= val - 1 using indicator approx)
        row = np.zeros(len(dirs))
        for i in range(len(dirs)):
            if x[i] > 0:
                row[i] = 1.0
        nogoods.append(LinearConstraint(row.reshape(1, -1), -np.inf, float(sum(1 for i in range(len(dirs)) if x[i] > 0) - 1)))
    return best
