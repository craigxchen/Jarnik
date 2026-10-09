"""NUMERICAL verification of a lift7 solution: build the 7 points (0,1,inf,a,b,c,d) as rational
functions of tau, find all collision fibres and count boundary contacts (class beta <=> each of the
56 splits exactly once)."""
import numpy as np, itertools, sys
from collections import Counter
from lift7b import Sys
from gen6 import mover_polys, MOV

def compose_y(P, yp, ym):
    """P(y) with y = (yp tau^2 - ym)/(tau^2 - 1): returns numerator in tau times (tau^2-1)^deg."""
    deg = len(P) - 1
    num = np.array([yp, 0, -ym], complex); den = np.array([1, 0, -1], complex)
    out = np.zeros(1, complex)
    for k, co in enumerate(P):   # co * y^(deg-k)
        term = np.array([1.0+0j])
        for _ in range(deg-k): term = np.polymul(term, num)
        for _ in range(k): term = np.polymul(term, den)
        out = np.polyadd(out, co*term)
    return out

def build(u, L):
    v0 = np.load('gen6_seed.npy')
    S = Sys(v0, np.random.default_rng(0)); S.L = L
    v, yp, ym, a = S.split(u)
    mp = mover_polys(v)
    pts = {'0': (np.array([0j]), np.array([1+0j])), '1': (np.array([1+0j]), np.array([1+0j])),
           'i': (np.array([1+0j]), np.array([0j]))}
    for m in MOV:
        k, q = mp[m]
        pts[m] = (k*compose_y(q['0'], yp, ym), compose_y(q['i'], yp, ym))
    pts['d'] = (a[:9], a[9:])
    return pts, (v, yp, ym, a)

def analyse(pts, tol=1e-6):
    names = list(pts)
    cand = []
    for p, q in itertools.combinations(names, 2):
        N = np.polysub(np.polymul(pts[p][0], pts[q][1]), np.polymul(pts[q][0], pts[p][1]))
        N = np.atleast_1d(N)
        while len(N) > 1 and abs(N[0]) < 1e-9*np.max(np.abs(N)): N = N[1:]
        if len(N) > 1: cand += list(np.roots(N))
    reps = []
    for x in cand:
        if not any(abs(x-r) < tol*max(1, abs(r)) for r in reps): reps.append(x)
    cnt = Counter(); bad = []
    P = set(names)
    for t0 in reps:
        vals = {p: (np.polyval(pts[p][0], t0), np.polyval(pts[p][1], t0)) for p in names}
        groups = []
        for p in names:
            for g in groups:
                n1, d1 = vals[p]; n2, d2 = vals[g[0]]
                sc = max(abs(n1), abs(d1))*max(abs(n2), abs(d2))
                if abs(n1*d2-n2*d1) < 1e-5*sc: g.append(p); break
            else: groups.append([p])
        big = [g for g in groups if len(g) >= 2]
        if len(big) != 1 or len(big[0]) > 5:
            bad.append((t0, big)); continue
        Sx = tuple(sorted(big[0])); Sc = tuple(sorted(P - set(Sx)))
        cnt[min(Sx, Sc)] += 1
    return len(reps), cnt, bad

if __name__ == '__main__':
    u = np.load(sys.argv[1]); L = np.load(sys.argv[2])
    pts, (v, yp, ym, a) = build(u, L)
    nrep, cnt, bad = analyse(pts)
    print('special fibres', nrep, 'bad fibres', len(bad), 'distinct splits', len(cnt),
          'max mult', max(cnt.values()) if cnt else None)
    print('branch points', yp, ym, ' |a| even-part check:', np.linalg.norm(a[1::2]), np.linalg.norm(a[0::2]))
    for t0, big in bad[:10]: print('   bad', t0, big)
