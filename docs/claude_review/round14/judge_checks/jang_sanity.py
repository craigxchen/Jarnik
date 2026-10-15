"""Sanity for jang.py: genuine level-1 labels (jglib.genuine) at all non-G moduli of each witness must give
PU-good and ANG-good systems, and their label-collisions must coincide with q | D_s(Z_v) (unit convention check)."""
import sys, os, glob, pickle
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jglib import two_sq, genuine, Zv, Ds, VK, vq
from jang import box_meets_circle
D = sys.argv[1]; bad = 0; checks = 0
for f in sorted(glob.glob(os.path.join(D, 'witness*.pkl'))):
    w = pickle.load(open(f, 'rb')); ps = w['ps']; pis = [two_sq(p) for p in ps]
    gl = {q: genuine(ps, q) for q in w['labels']}
    for v in VK(len(ps), w['K']):
        z = Zv(pis, v); N = z[0]**2 + z[1]**2; ms = []
        for s in range(4):
            d = Ds(z, s); mm = 1
            for q in w['G']: mm *= q**vq(d, q)
            for q, (m, L) in gl.items():
                coll = sum(a*b for a, b in zip(v, L)) % m == s*(m//4)
                checks += 1
                if coll != (d % q == 0): bad += 1
                if coll: mm *= q
            ms.append(mm)
        m0, m1, m2, m3 = ms
        if not box_meets_circle(Fr(2*m2*m2, N)-1, 1-Fr(2*m0*m0, N), Fr(m3*m3, N)-1, 1-Fr(m1*m1, N)): bad += 1
print('checks', checks, 'mismatches/ANG failures', bad)
