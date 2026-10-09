"""NUMERICAL: verify the k=6 Klein-ansatz solutions (klein6.py / klein6t.py): collision fibres of the
6-point configuration (0,1,inf,a,b,c) as functions of t; each of the 25 splits must occur exactly once."""
import numpy as np, itertools
from collections import Counter
from klein6 import movers_t, collisions, build

def analyse(M):
    names, ev = collisions(M)
    cnt = Counter(); bad = 0; P = set(names)
    for t0, big in ev:
        if len(big) != 1 or len(big[0]) > len(P)-2: bad += 1; continue
        S = tuple(sorted(big[0])); Sc = tuple(sorted(P - set(S)))
        cnt[min(S, Sc)] += 1
    return len(ev), bad, cnt

u = np.load('klein6_sol_3.npy'); zc = np.load('klein6_zc.npy')
nf, bad, cnt = analyse(movers_t(u, zc))
print('klein6 seed 3: fibres', nf, 'bad', bad, 'distinct splits', len(cnt), 'max mult', max(cnt.values()))
# real/complex nature
v1, y1, v2, y2, kA, kB, v0, y0 = build(u, zc)
print('  max |Im| of parameters:', np.max(np.abs(np.imag(u))))
# klein6t solutions: convert to the klein6 representation is not needed; analyse via t-formulation
sols = np.load('klein6t_sols.npy'); fixed = np.load('klein6t_fixed.npy')
good = 0
for r in sols:
    t = np.zeros((3, 2), complex); t[:, 0] = fixed; t[:, 1] = r[0:3]
    kA, kB = r[4], r[5]
    v = t + 1/t; y = t**2; z = y + 1/y
    PA = [np.poly(v[x]) for x in range(3)]; PB = [np.poly(y[x]) for x in range(3)]
    def sub_v(Pp):
        deg = len(Pp)-1; out = np.zeros(1, complex)
        for k, co in enumerate(Pp):
            term = np.array([1.0+0j])
            for _ in range(deg-k): term = np.polymul(term, [1, 0, 1])
            term = np.polymul(term, np.concatenate([[1.0], np.zeros(k)]))
            out = np.polyadd(out, co*term)
        return out
    def sub_y(Pp):
        out = np.zeros(1, complex)
        for co in Pp: out = np.polyadd(np.polymul(out, [1, 0, 0]), [co])
        return out
    zc = z[:, 0]
    zpoly = lambda rr: np.array([1, 0, -rr, 0, 1], complex)
    M = {'a': (kA*sub_v(PA[0]), sub_v(PA[2])), 'b': (kB*sub_y(PB[0]), sub_y(PB[2])),
         'c': ((zc[1]-zc[2])*zpoly(zc[0]), (zc[1]-zc[0])*zpoly(zc[2]))}
    nf, bad, cnt = analyse(M)
    ok = (bad == 0 and len(cnt) == 25 and max(cnt.values()) == 1)
    good += ok
print('klein6t solutions:', len(sols), 'nondegenerate (25 splits once each):', good)
