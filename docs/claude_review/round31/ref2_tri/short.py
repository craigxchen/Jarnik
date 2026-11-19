"""ref2_tri: Astra's computed triangle as a genuine short triple; and the 417-digit 4-point cluster
(construction recipe copied from astra_free_II/near_square_check.py, analysed with own check_triple)."""
import sys; sys.path.insert(0, '.')
from genuine import *
N = 131490079981
pts = points_on(N)
labels = {6451601, 28879877, 92794693}
found = set()
for i in range(len(pts)):
    for j in range(len(pts)):
        if i == j: continue
        A = half(pts[i], pts[j], N)
        if A and nm(A) in labels and abs(A[1]) <= 2 and abs(A[0]) > 1000: found.add(pts[i])
found = sorted(found, key=lambda z: math.atan2(z[1], z[0]))
# find 3 mutually linked points
res = []
for t in combinations(found, 3):
    ok = all((lambda A: A and nm(A) in labels and abs(A[1]) <= 2)(half(p, q, N)) for p, q in combinations(t, 2))
    if ok:
        for rot in [(0,1,2),(1,2,0),(2,0,1)]:
            r = check_triple(t[rot[0]], t[rot[1]], t[rot[2]], N)
            res.append(r)
        break
print('short genuine triangle at N=131490079981, three anchors:', res)
# Vieta reflection in Astra's normalisation
a,b,c,u,v,w = 1417,20381,4553,2,1,1
A_,B_,C_ = a*u*u, b*v*v, c*w*w
U = (u*v*w)**2
print('A,B,C,U', A_,B_,C_,U, 'Heron', (A_+B_-C_)**2 == 4*A_*B_-4*U, 'abc==N', a*b*c==N)
Bp = 2*(A_+C_) - B_
print('reflected', Bp, 'N mod', N % Bp, 'check product', Bp*B_ == (A_-C_)**2+4*U)

# 417-digit cluster
from fractions import Fraction
n = 10**30+1
x, y, z, uu = -3, -2, 6, -6
blocks = []
for p_, q_, r_ in [(x, y, z), (x, z, y), (y, z, x)]:
    blocks.append({(0, 1), (uu, p_-1), (uu, q_-1), (uu*(1+r_)//r_, r_)})
roots = sorted(set.union(*blocks))
P = []
for chosen in [set(), *blocks]:
    f = [(n-4*re, (-4 if (re, im) in chosen else 4)*im) for re, im in roots]
    P.append(reduce(mul, f, (1, 0)))
g = reduce(ggcd, P)
P = [ (lambda q: q)( (lambda zz: (zz[0]//nm(g), zz[1]//nm(g)))(mul(p, cj(g))) ) for p in P]
NN = nm(P[0]); assert all(nm(p) == NN for p in P)
before = stats['ok']; out = []
for t in combinations(range(4), 3):
    for rot in [(0,1,2),(1,2,0),(2,0,1)]:
        r = check_triple(P[t[rot[0]]], P[t[rot[1]]], P[t[rot[2]]], NN)
        if r: out.append((r[3], r[4], r[5], r[7]))
print('417-digit cluster: digits', len(str(NN)), 'anchored triples ok', stats['ok']-before, 'of 12; (u,v,w,eps):', out[:4])
print(stats)
