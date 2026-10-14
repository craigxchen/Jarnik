"""Which planes inside the Fano 3-space U_(1,2,3) contribute degree-6 elements of J(A) outside J_w ?"""
import itertools, math, pickle
import numpy as np
from m8lib import *
A, cols, H, JA, W = pickle.load(open('m8_step1.pkl','rb'))
gw = pickle.load(open('m8_gw.pkl','rb'))
rows_w6 = ideal_degree_rows(gw, 6); base = dim_exact(rows_w6, 6)
ca, cb, cc = cols[1], cols[2], cols[3]
z = tuple((x+y+w)//2 for x,y,w in zip(ca,cb,cc))
basis = (ca, cb, z)
R = 3
cache = {}
normals = [v for v in itertools.product(range(-R, R+1), repeat=3)
           if any(v) and [t for t in v if t][0] > 0 and math.gcd(*[abs(t) for t in v if t]) == 1]
def plane_basis(nv):
    cands = [v for v in itertools.product(range(-R-1, R+2), repeat=3) if any(v) and sum(p*q for p,q in zip(v,nv)) == 0]
    cands.sort(key=lambda v: sum(abs(t) for t in v))
    b1 = cands[0]
    for v in cands[1:]:
        cr = (b1[1]*v[2]-b1[2]*v[1], b1[2]*v[0]-b1[0]*v[2], b1[0]*v[1]-b1[1]*v[0])
        if any(cr): return b1, v
Aa = np.array(A)
allg = []
for nv in normals:
    b1, b2 = plane_basis(nv)
    F = tuple(tuple(sum(bb[i]*basis[i][x] for i in range(3)) for x in range(7)) for bb in (b1, b2))
    g, rho, ny = proj_ideal_gens(A, F, 8)
    if 6 in g:
        pb = [pullback(f, F, cache) for f in g[6]]
        gain = dim_exact(rows_w6 + pb, 6) - base
        Y = sorted(set(map(tuple, (Aa @ np.array(F).T).tolist())))
        print(f"normal {nv} (coords wrt c_a,c_b,z) basis {b1},{b2}: F={F}  |Y|={ny} rho(Y)={rho} #deg6 gens={len(g[6])} gain over J_w={gain}")
        allg += pb
print("total gain of all these planes over J_w in degree 6:", dim_exact(rows_w6 + allg, 6) - base)
