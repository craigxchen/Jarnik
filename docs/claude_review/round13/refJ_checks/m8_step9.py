"""Admissibility (rho(pi_F(A)) < m = 8) of the contributing families, and structure of the images.
(1) column pair image = full grid [-2,2]^2 ; non-line triple image = full cube [-2,2]^3 (so J = (y_i^5), pullbacks are width powers).
(2) Fano-line triple image: |Y| = 57, contained in the parity lattice y_a+y_b+y_{a+b} even.
(3) planes P in the Fano 3-space U through z (normals (n1,n2,0) wrt basis c_a,c_b,z): rho(Y_P) to degree 10, and the
    degree-6 gain over J_w of the ADMISSIBLE ones (rho < 8), all 7 lines; and of all of them (coeff bound R)."""
import itertools, math, pickle, time
import numpy as np
from m8lib import *
A, cols, H, JA, W = pickle.load(open('m8_step1.pkl','rb'))
gw = pickle.load(open('m8_gw.pkl','rb'))
Aa = np.array(A)
rows_w6 = ideal_degree_rows(gw, 6); base = dim_exact(rows_w6, 6)
img = lambda F: sorted(set(map(tuple, (Aa @ np.array(F).T).tolist())))
Yp = img((cols[1], cols[2])); print("column pair (1,2): |Y| =", len(Yp), " full grid:", set(Yp) == set(itertools.product(range(-2,3), repeat=2)))
print("all 21 column pairs full grid:", all(set(img((cols[a],cols[b]))) == set(itertools.product(range(-2,3), repeat=2)) for a,b in itertools.combinations(range(1,8),2)))
lines = [(1,2,3),(1,4,5),(1,6,7),(2,4,6),(2,5,7),(3,4,7),(3,5,6)]
nonl = [t for t in itertools.combinations(range(1,8),3) if t not in lines]
print("all 28 non-line triples full cube:", all(set(img(tuple(cols[a] for a in t))) == set(itertools.product(range(-2,3), repeat=3)) for t in nonl))
Yf = img(tuple(cols[a] for a in (1,2,3)))
par = [y for y in itertools.product(range(-2,3), repeat=3) if sum(y) % 2 == 0]
print("Fano triple (1,2,3): |Y| =", len(Yf), " subset of parity lattice box (63 pts):", set(Yf) <= set(par), " missing:", sorted(set(par) - set(Yf)))
R = 4; cache = {}
adm_all = []; all_all = []
for (a,b,c) in lines:
    ca, cb, cc = cols[a], cols[b], cols[c]
    z = tuple((x+y+w)//2 for x,y,w in zip(ca,cb,cc))
    for n1, n2 in [(n1,n2) for n1 in range(0,R+1) for n2 in range(-R,R+1) if (n1,n2)!=(0,0) and math.gcd(n1,abs(n2))==1 and (n1>0 or n2>0)]:
        # plane {alpha c_a + beta c_b + gamma z : n1 alpha + n2 beta = 0} = span(z, n2 c_a - n1 c_b)
        v = tuple(n2*p - n1*q for p, q in zip(ca, cb))
        F = (z, v)
        g, rho, ny = proj_ideal_gens(A, F, 10)
        pb = [pullback(f, F, cache) for f in g.get(6, [])]
        all_all += pb
        if rho is not None and rho < 8: adm_all += pb
        if (a,b,c) == (1,2,3): print(f"  line (1,2,3) plane span(z, {n2}c_a-{n1}c_b): |Y|={ny} rho={rho} #deg6 gens={len(g.get(6,[]))}")
print("degree-6 gain over J_w, all 7 lines: admissible pencil planes (rho<8):", dim_exact(rows_w6 + adm_all, 6) - base,
      "; all pencil planes |n|<=%d:" % R, dim_exact(rows_w6 + all_all, 6) - base)
