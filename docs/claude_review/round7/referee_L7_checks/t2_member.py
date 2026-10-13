# Theorem B class of the member through x, and comp for sigma=id (should be 2 = lambda(gamma_A)).
import sys, time; sys.path.insert(0, '.')
from reflib import *
cube = {(i,j,k): meet3(i,j,k) for i in (1,2) for j in (3,4) for k in (5,6)}
Ekeys = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
webm = web([cube[k] for k in Ekeys])
def show(cls):
    out = []
    for S, d in sorted(cls.items(), key=lambda kv: (len(kv[0]), sorted(kv[0]))):
        out.append(('-'.join(LAB[i] for i in sorted(S)), d))
    return out
for x in [(6, -6, -3, -5, 2, 0, 1), (3, -2, 4, 0, 2, -5, -6), (1, 7, -4, 2, -9, 5, 3)]:
    V, data = member_through(x, webm)
    # check: at s=0 the configuration is projectively x
    v0 = [(int(V[a][0].coeffs()[0]) if V[a][0] != 0 else 0, int(V[a][1].coeffs()[0]) if V[a][1] != 0 else 0) for a in range(7)]
    from fractions import Fraction as Fr
    def cr(z, a, b, c, d):
        D = lambda u, v: u[0]*v[1]-u[1]*v[0]
        return Fr(D(z[a],z[c])*D(z[b],z[d]), D(z[a],z[d])*D(z[b],z[c]))
    xz = [(Fr(t), 1) for t in x]
    xz = [(Fr(t).numerator, Fr(t).denominator) for t in x]
    same = all(cr(v0, *k) == cr(xz, *k) for k in FOURS)
    cls, good, info = class_and_goodness(V)
    print('x =', x, ' member passes through x (cross-ratios at s=0 agree):', same)
    print('  good:', good, info)
    print('  class:', show(cls))
    pr = sum(d for S, d in cls.items() if len(S) == 2); tr = sum(d for S, d in cls.items() if len(S) == 3)
    print('  pair', pr, 'triple', tr, '(B+8K)=', tr - Fr(5, 3) * pr)
    t0 = time.time()
    Br, dBr, _ = FA_brackets(V, webm)
    print('  comp(F_A, member) =', comp_count(Br, dBr), ' time %.1fs' % (time.time() - t0))
