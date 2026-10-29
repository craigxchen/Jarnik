"""Zonotope orders on slab boxes (own code).
zeta_pair: exhaustive DFS over L_pair generators (with multiplicity).  zeta_all: exact LP-dual upper bound over all
primitive v in {-2..2}^d, then exhaustive DFS over the generators whose dual cost leaves room for the target order.
usage: zono_v.py SPEC [mode]   mode = pair | all (default all = both)
"""
import sys, json, time, itertools, math
from fractions import Fraction
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks')
from vlib import *

def functionals(n, pats, A):
    """rows: integer functionals on k' (coordinates, base coordinate sigma, slabs); widths over A."""
    d = n - 1
    F = []
    for x in range(d):
        v = [0] * d; v[x] = 1; F.append(v)
    F.append([1] * d)  # sigma = -k_base
    for S in pats:
        # kfull(S) = sum_{x in S, x<d} k_x - [d in S] * sum_all k_x
        v = [0] * d
        for x in range(d):
            v[x] = (1 if x in S else 0) - (1 if d in S else 0)
        F.append(v)
    F = np.array(F, dtype=np.int64)
    vals = A.astype(np.int64) @ F.T
    W = vals.max(axis=0) - vals.min(axis=0)
    return F, W

def candidates(d, F, W, pair_only=False):
    if pair_only:
        return np.array(lpair_classes(d), dtype=np.int64)
    C = np.array(list(itertools.product(range(-2, 3), repeat=d)), dtype=np.int64)
    # sign-normalise: first nonzero positive
    nz = (C != 0)
    first = np.argmax(nz, axis=1)
    sgn = C[np.arange(len(C)), first]
    C = C[(sgn > 0)]
    g = np.zeros(len(C), dtype=np.int64)
    for x in range(d): g = np.gcd(g, C[:, x])
    C = C[g == 1]
    cost = np.abs(C @ F.T)
    C = C[(cost <= W[None, :]).all(axis=1)]
    return C

def lp_dual(F, W, C):
    """min W.y  s.t.  |C F^T| y >= 1, y >= 0.  Returns exact Fractions y (feasible, verified) and bound B."""
    from scipy.optimize import linprog
    Aub = -np.abs(C @ F.T).astype(float)
    res = linprog(W.astype(float), A_ub=Aub, b_ub=-np.ones(len(C)), bounds=[(0, None)] * len(W), method='highs')
    assert res.status == 0, res.message
    y = [Fraction(v).limit_denominator(2000) if v > 1e-12 else Fraction(0) for v in res.x]
    D = math.lcm(*[f.denominator for f in y])
    Y = np.array([int(f * D) for f in y], dtype=np.int64)
    cv = np.abs(C @ F.T) @ Y           # = D * c_v, exact integers
    m = int(cv.min())                    # rescale so that min c_v = 1 exactly
    # exact: c_v = cv / m >= 1 for all candidates; bound B = W.Y / m
    B = Fraction(int(W @ Y), m)
    return Y, m, B, cv

class DFS:
    def __init__(self, A, gens, F, Y, m, costs_excess=None):
        self.A = A; self.N = len(A)
        idx = {tuple(k.tolist()): i for i, k in enumerate(A)}
        self.gens = gens
        self.nbr = []
        for v in gens:
            nb = np.full(self.N, self.N, dtype=np.int64)
            for i, k in enumerate(A):
                t = tuple((k.astype(np.int64) + v).tolist())
                if t in idx: nb[i] = idx[t]
            self.nbr.append(nb)
        self.vals = A.astype(np.int64) @ F.T  # N x nF
        self.Y = Y; self.m = m
        self.excess = costs_excess  # per generator: D*c_v - m (>= 0), or None
        self.nodes = 0

    def bound(self, K):
        v = self.vals[K]
        wid = v.max(axis=0) - v.min(axis=0)
        return int(wid @ self.Y) // self.m   # sum_w c_w <= W(K).y, c_w >= 1

    def run(self, target=None, budget=None):
        """max order (or first zonotope of order >= target). budget: max total excess (in units of 1/m... times D)."""
        self.best = -1; self.wit = None; self.target = target
        K = np.ones(self.N, dtype=bool)
        self._dfs(K, 0, [], 0, budget)
        return self.best, self.wit

    def _dfs(self, K, start, chosen, exc, budget):
        self.nodes += 1
        depth = len(chosen)
        if depth > self.best:
            self.best = depth; self.wit = (int(np.argmax(K)), list(chosen))
        if self.target is not None and self.best >= self.target: return True
        for i in range(start, len(self.gens)):
            e2 = exc
            if budget is not None:
                e2 = exc + int(self.excess[i])
                if e2 > budget: continue
            Kx = np.append(K, False)
            K2 = K & Kx[self.nbr[i]]
            if not K2.any(): continue
            need = (self.target if self.target is not None else self.best + 1) - depth - 1
            if need > 0 and self.bound(K2) < need: continue
            chosen.append(i)
            if self._dfs(K2, i, chosen, e2, budget): return True
            chosen.pop()
        return False

def verify_witness(A, k0, vs):
    """exact: all 2^t subset sums in A; moments of order < t vanish on std monomials (enough on the cube); some order-t
    moment nonzero. Returns (ok, |c|_inf, support size)."""
    S = set(tuple(k.tolist()) for k in A)
    t = len(vs)
    coef = {}
    for T in itertools.product([0, 1], repeat=t):
        pt = tuple(int(x) for x in (np.array(k0) + sum((np.array(v) * b for v, b in zip(vs, T)), np.zeros(len(k0), dtype=np.int64))))
        if pt not in S: return (False, 'point outside A', pt)
        coef[pt] = coef.get(pt, 0) + (-1) ** (t - sum(T))
    coef = {k: c for k, c in coef.items() if c}
    d = len(k0)
    by = std_monos(d, t)
    pts = np.array(list(coef.keys()), dtype=np.int8); cs = np.array(list(coef.values()), dtype=np.int64)
    for e in range(t):
        Em = evalmat(pts, by[e]).astype(np.int64)
        mom = cs @ Em
        if (mom != 0).any(): return (False, 'moment nonzero at degree', e)
    # order-t: <k,u>^t moment: sum c_k (k.u)^t at random integer u must equal t! prod (v.u) != 0
    u = np.array([3 ** x + 7 for x in range(d)], dtype=object)
    s = sum(c * int(np.dot(np.array(k, dtype=object), u)) ** t for k, c in coef.items())
    pr = math.factorial(t)
    for v in vs: pr *= int(np.dot(np.array(v, dtype=object), u))
    return (s == pr and s != 0, max(abs(c) for c in coef.values()), len(coef))

if __name__ == '__main__':
    spec = sys.argv[1]; mode = sys.argv[2] if len(sys.argv) > 2 else 'both'
    t0 = time.time()
    n, pats, L, A = spec_box(spec); d = n - 1
    F, W = functionals(n, pats, A)
    out = {'spec': spec, 'size': len(A)}
    # pair
    Cp = candidates(d, F, W, pair_only=True)
    Yp, mp, Bp, cvp = lp_dual(F, W, Cp)
    dp = DFS(A, Cp, F, Yp, mp)
    zp, wp = dp.run()
    ok = verify_witness(A, A[wp[0]].astype(int).tolist(), [Cp[i].tolist() for i in wp[1]])
    out.update(zeta_pair=zp, lp_pair=str(Bp), nodes_pair=dp.nodes, wit_pair=[A[wp[0]].tolist(), [Cp[i].tolist() for i in wp[1]]], wit_pair_ok=str(ok))
    log(spec, 'zeta_pair', zp, 'LP_pair', Bp, 'nodes', dp.nodes, 'witness ok', ok)
    if mode != 'pair':
        C = candidates(d, F, W)
        Y, m, B, cv = lp_dual(F, W, C)
        T = int(B)  # floor
        out.update(ncand=len(C), lp_all=str(B))
        log('candidates', len(C), 'LP_all bound', B, '=', float(B))
        if T <= zp:
            out['zeta_all'] = zp; out['zeta_all_reason'] = 'floor(LP_all) = %d <= zeta_pair' % T
        else:
            # search orders T, T-1, ..., zp+1 over generators with c_v - 1 <= B - t
            za = zp
            for t in range(T, zp, -1):
                slack = B - t  # exact Fraction; excess_v = c_v - 1 = (cv - m)/m <= slack  <=> cv - m <= slack*m
                exc = cv - m
                lim = slack * m
                sel = [i for i in range(len(C)) if exc[i] <= lim]
                G = C[sel]
                # budget in integer units: sum (cv - m) <= slack*m  ->  multiply by denominator of slack*m
                den = Fraction(lim).denominator
                excI = (exc[sel] * den).astype(np.int64)
                budget = int(lim * den)
                dd = DFS(A, G, F, Y, m, costs_excess=excI)
                b, w = dd.run(target=t, budget=budget)
                log('order', t, 'generators allowed', len(G), 'nodes', dd.nodes, 'found' if b >= t else 'none')
                out['search_%d' % t] = {'ngen': len(G), 'nodes': dd.nodes, 'found': b >= t,
                                        'gens_are_pair': all(tuple(g) in set(lpair_classes(d)) for g in G.tolist())}
                if b >= t:
                    za = t
                    out['wit_all'] = [A[w[0]].tolist(), [G[i].tolist() for i in w[1]]]
                    out['wit_all_ok'] = str(verify_witness(A, A[w[0]].astype(int).tolist(), [G[i].tolist() for i in w[1]]))
                    break
            out['zeta_all'] = za
        log(spec, 'zeta_all', out['zeta_all'])
    out['time'] = round(time.time() - t0, 1)
    with open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks/zeta.jsonl', 'a') as f:
        f.write(json.dumps(out) + '\n')
