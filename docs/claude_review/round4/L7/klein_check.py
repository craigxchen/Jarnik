"""NUMERICAL verification of a Klein-ansatz solution: build the 7 points as rational functions
of t, find every collision fibre, and list the boundary divisors met (with multiplicity)."""
import numpy as np, sys, itertools
from klein import unpack, maps

def ratfun_t(u):
    v, y, z, PA, PB, QC, QD, kA, kB, kC, kD, vcd, ycd = maps(u)
    # a = kA PA0(v)/PAinf(v), v = t + 1/t -> multiply num/den by t^4
    def sub_v(P):  # P(t+1/t)*t^deg
        deg = len(P)-1; out = np.zeros(1)
        for k, co in enumerate(P):
            # co * (t+1/t)^(deg-k) * t^deg = co*(t^2+1)^(deg-k) t^k
            term = np.array([1.0+0j])
            for _ in range(deg-k): term = np.polymul(term, [1, 0, 1])
            term = np.polymul(term, np.concatenate([[1.0], np.zeros(k)]))
            out = np.polyadd(out, co*term)
        return out
    def sub_y(P):  # P(t^2)
        out = np.zeros(1)
        for co in P:
            out = np.polyadd(np.polymul(out, [1, 0, 0]), [co])
        return out
    def sub_z(Q):  # Q(t^2+t^-2)*t^(2deg) = co (t^4+1)^(deg-k) t^(2k)
        deg = len(Q)-1; out = np.zeros(1)
        for k, co in enumerate(Q):
            term = np.array([1.0+0j])
            for _ in range(deg-k): term = np.polymul(term, [1, 0, 0, 0, 1])
            term = np.polymul(term, np.concatenate([[1.0], np.zeros(2*k)]))
            out = np.polyadd(out, co*term)
        return out
    a = (kA*sub_v(PA[0]), sub_v(PA[2]))
    b = (kB*sub_y(PB[0]), sub_y(PB[2]))
    c = (kC*sub_z(QC[0]), sub_z(QC[2]))
    d = (kD*sub_z(QD[0]), sub_z(QD[2]))
    return {'a': a, 'b': b, 'c': c, 'd': d}

def analyse(u, tol=1e-6, verbose=True):
    M = ratfun_t(u)
    pts = {'0': (np.array([0.0+0j]), np.array([1.0+0j])), '1': (np.array([1.0+0j]), np.array([1.0+0j])),
           'i': (np.array([1.0+0j]), np.array([0.0+0j]))}
    pts.update(M)
    names = list(pts)
    # collision t-values: roots of N1*D2 - N2*D1 for each pair
    cand = []
    for p, q in itertools.combinations(names, 2):
        N = np.polysub(np.polymul(pts[p][0], pts[q][1]), np.polymul(pts[q][0], pts[p][1]))
        N = np.trim_zeros(N, 'f')
        # strip tiny leading coefficients
        while len(N) > 1 and abs(N[0]) < 1e-9*np.max(np.abs(N)): N = N[1:]
        cand += list(np.roots(N))
    # also t = 0 and t = infinity: check separately
    cand = np.array(cand)
    # cluster candidate t-values
    reps = []
    for x in cand:
        if not any(abs(x - r) < tol*max(1, abs(r)) for r in reps):
            reps.append(x)
    events = []
    for t0 in reps:
        vals = {}
        for p in names:
            Nn, Dd = pts[p]
            nv = np.polyval(Nn, t0); dv = np.polyval(Dd, t0)
            vals[p] = (nv, dv)
        # group by projective equality
        groups = []
        for p in names:
            placed = False
            for g in groups:
                q = g[0]
                n1, d1 = vals[p]; n2, d2 = vals[q]
                sc = max(abs(n1), abs(d1))*max(abs(n2), abs(d2))
                if abs(n1*d2 - n2*d1) < 1e-5*sc:
                    g.append(p); placed = True; break
            if not placed: groups.append([p])
        big = [g for g in groups if len(g) >= 2]
        events.append((t0, big))
    return events

if __name__ == '__main__':
    u = np.load(sys.argv[1])
    ev = analyse(u)
    from collections import Counter
    cnt = Counter()
    bad = 0
    for t0, big in ev:
        if len(big) != 1:
            print('fibre with clusters', big, 't=', t0); bad += 1; continue
        S = frozenset(big[0])
        if len(S) >= 6:
            print('fibre with', len(S), '-cluster', sorted(S)); bad += 1; continue
        Sc = frozenset(['0', '1', 'i', 'a', 'b', 'c', 'd']) - S
        key = min(tuple(sorted(S)), tuple(sorted(Sc)))
        cnt[key] += 1
    print('number of special fibres', len(ev), 'bad', bad)
    print('distinct splits met', len(cnt), 'max multiplicity', max(cnt.values()) if cnt else None)
    allsplits = set()
    P7 = ['0', '1', 'i', 'a', 'b', 'c', 'd']
    for r in [2, 3]:
        for S in itertools.combinations(P7, r):
            Sc = tuple(sorted(set(P7)-set(S)))
            allsplits.add(min(tuple(sorted(S)), Sc))
    print('total splits', len(allsplits), 'missing', [s for s in allsplits if s not in cnt][:10])
