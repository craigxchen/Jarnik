import sys, random, itertools
import sympy as sp
from sympy import Rational as R
from elam import deg_EL, basis_quadrics, t
from star import point_of, LAB
from arr import PLANES, CHOSEN
pts = [point_of([PLANES[i], PLANES[j], PLANES[k]]) for (i,j,k) in CHOSEN]
qb = basis_quadrics(pts); assert len(qb) == 4
def rnd(): 
    while True:
        n = random.randint(-60, 60)
        if n: return R(n, random.randint(1, 25))
def good_curve(T):
    while True:
        c = [[rnd() for _ in range(4)] for a in range(7)]
        if T:
            a0 = T[0]; z0 = c[a0][1]/c[a0][3]
            for b in T[1:]: c[b][1] = z0*c[b][3]
        if any(ci[0]*ci[3] - ci[1]*ci[2] == 0 for ci in c): continue
        tt = R(7, 13); vals = [(ci[0]*tt + ci[1])/(ci[2]*tt + ci[3]) for ci in c]
        if len(set(vals)) < 7: continue
        return [(ci[0]*t + ci[1], ci[2]*t + ci[3]) for ci in c]
T = tuple(int(s) for s in sys.argv[1].split(',')) if len(sys.argv) > 1 and sys.argv[1] != 'mu' else ()
random.seed(hash(T) % 100000 + 17)
out = []
for rep in range(2):
    d, _ = deg_EL(PLANES, qb, good_curve(T)); out.append(d)
print('T', [LAB[a] for a in T] if T else 'mu', 'E_Lambda degrees', out, flush=True)
