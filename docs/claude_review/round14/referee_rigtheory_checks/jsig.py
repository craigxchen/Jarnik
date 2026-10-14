"""Exact signature sets for the joint problem (referee's own). Signature at q: frozenset of (row, a) meaning
(v_row, s_row) collides at levels 1..a (we store max a per row). Labels outside the eps-class at q0 flagged."""
import sys, math, itertools, pickle, numpy as np
from glib import *
from forced import forced_set
def canon(v):
    for x in v:
        if x: return x > 0
def U1_orbit_reps(m, par):
    seen = set(); reps = []
    U = [u for u in range(1, m, 4) if math.gcd(u, m) == 1]
    for x in range(par, m, 2):
        if x in seen: continue
        reps.append(x)
        for u in U: seen.add((u*x) % m)
    return reps
def sigs(ps, K, q, lev=1, eps_class=True):
    r = len(ps); pis = [two_sq(p) for p in ps]; V = [v for v in VK(r, K) if canon(v)]
    nz = [Zv(pis, v)[0]**2 + Zv(pis, v)[1]**2 for v in V]
    m1 = q-1 if q % 4 == 1 else q+1; n = m1//4; mt = m1*q**(lev-1)
    ml = [m1*q**(a-1) for a in range(1, lev+1)]
    # affordable depth per (v,s): q^(2a) <= nz*2^[s odd]
    dep = np.zeros((len(V), 4), np.int64)
    for i in range(len(V)):
        for s in range(4):
            d = 0
            while d < lev and q**(2*(d+1)) <= nz[i]*(2 if s % 2 else 1): d += 1
            dep[i, s] = d
    nuv = nu(ps, q); Va = np.array(V, np.int64)
    m0, L = genuine(ps, q)
    cls = set()
    if n > 1:
        for eps in itertools.product([1, -1], repeat=r):
            for u in range(1, n):
                if math.gcd(u, n) == 1: cls.add(tuple((u*e*x) % n for e, x in zip(eps, L)))
    out = {}   # signature -> flags (1 class, 2 nonclass)
    ch = [np.arange(nuv[j], mt, 2) for j in range(r)]
    rest = np.array(np.meshgrid(*ch[1:], indexing='ij')).reshape(r-1, -1).T
    for a0 in U1_orbit_reps(mt, nuv[0]):
        grid = np.hstack([np.full((len(rest), 1), a0), rest])
        t = (grid @ Va.T) % mt
        t1 = t % m1
        coll = (t1 % n == 0)
        s = np.where(coll, (t1//n) % 4, 0)
        depth = np.zeros_like(t)
        for a in range(lev):
            depth += ((t % ml[a]) % (ml[a]//4) == 0) & (depth == a)
        depth = np.where(coll, depth, 0)
        aff = (depth <= dep[np.arange(len(V))[None, :], s]).all(axis=1)
        idxs = np.nonzero(aff)[0]
        for k in idxs:
            cv = np.nonzero(coll[k])[0]
            sig = tuple(sorted((int(i)*4 + int(s[k, i]), int(depth[k, i])) for i in cv))
            if n > 1:
                c = tuple(int(x) % n for x in grid[k])
                f = 1 if c in cls else 2
            else: f = 1
            out[sig] = out.get(sig, 0) | f
    return out
def minimal(S):
    S = sorted(S, key=len); res = []
    for s in S:
        ss = set(s)
        # dominated if some res r has r subset of s (with levels <=)
        dom = False
        for r_ in res:
            if all(any(row == rr and lv >= ll for (row, lv) in s) for (rr, ll) in r_): dom = True; break
        if not dom: res.append(s)
    return res
if __name__ == "__main__":
    ps = [int(x) for x in sys.argv[1].split(",")]; K = int(sys.argv[2]); levs = eval(sys.argv[3]); q0 = int(sys.argv[4])
    F, _ = forced_set(ps, K, 2000)
    data = {}
    for q in F:
        lev = levs.get(q, 1)
        o = sigs(ps, K, q, lev)
        allmin = minimal(list(o)); nonc = minimal([s for s, f in o.items() if f & 2])
        data[q] = dict(lev=lev, allmin=allmin, nonclass=nonc)
        print(q, "lev", lev, "sigs", len(o), "min", len(allmin), "min-nonclass", len(nonc), flush=True)
    pickle.dump((F, data), open(f"jsig_K{K}_{q0}.pkl", "wb"))
