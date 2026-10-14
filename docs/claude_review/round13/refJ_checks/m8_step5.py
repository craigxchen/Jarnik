"""Do PAIRS of general integer directions (2-dim rational subspaces U' inside the Fano 3-spaces U) fill the degree-6 gap?
J_{U'} := pullback of J(pi_{U'}(A)) depends only on the rational plane U' (a rational change of basis acts on the
image by an invertible linear map, and leading-form ideals are equivariant).  Planes in U = span(c_a,c_b,z) are
parametrised by primitive normals nvec in [-R,R]^3 w.r.t. the basis (c_a, c_b, z)."""
import time, pickle, itertools, math, sys
from m8lib import *
A, cols, H, JA, W = pickle.load(open('m8_step1.pkl','rb'))
gw = pickle.load(open('m8_gw.pkl','rb'))
R = int(sys.argv[1]) if len(sys.argv) > 1 else 3
DEG = 6
cache = {}
lines = [(1,2,3),(1,4,5),(1,6,7),(2,4,6),(2,5,7),(3,4,7),(3,5,6)]
def plane_basis(nv):
    # two independent integer vectors orthogonal to nv in Z^3
    cands = [v for v in itertools.product(range(-R-1, R+2), repeat=3) if any(v) and sum(p*q for p,q in zip(v,nv)) == 0]
    cands.sort(key=lambda v: sum(abs(t) for t in v))
    b1 = cands[0]
    for v in cands[1:]:
        cr = (b1[1]*v[2]-b1[2]*v[1], b1[2]*v[0]-b1[0]*v[2], b1[0]*v[1]-b1[1]*v[0])
        if any(cr): return b1, v
normals = [v for v in itertools.product(range(-R, R+1), repeat=3)
           if any(v) and [t for t in v if t][0] > 0 and math.gcd(*[abs(t) for t in v if t]) == 1]
t0 = time.time()
extra = {}
nplanes = 0; contributing = 0; rhos = {}
for (a,b,c) in lines:
    ca, cb, cc = cols[a], cols[b], cols[c]
    z = tuple((x+y+w)//2 for x,y,w in zip(ca,cb,cc))
    basis = (ca, cb, z)
    for nv in normals:
        b1, b2 = plane_basis(nv)
        F = tuple(tuple(sum(bb[i]*basis[i][x] for i in range(3)) for x in range(7)) for bb in (b1, b2))
        g, rho, ny = proj_ideal_gens(A, F, DEG)
        nplanes += 1
        rhos[rho] = rhos.get(rho, 0) + 1
        for e, fl in g.items():
            extra.setdefault(e, []).extend(pullback(f, F, cache) for f in fl)
print(f"R={R}: {nplanes} planes inside the 7 Fano 3-spaces; rho histogram {rhos}; generator counts {{e: len}} = { {e: len(l) for e,l in extra.items()} } [{time.time()-t0:.0f}s]", flush=True)
G = {}
for gg in (gw, extra):
    for e, l in gg.items(): G.setdefault(e, []).extend(l)
for i in range(DEG+1):
    rows = ideal_degree_rows(G, i)
    print(f"  deg {i}: dim_p = {dim_modp(rows, i)}  (J_w: {[0,0,0,7,49,203,616][i]}, J(A): {JA[i]})", flush=True)
rows = ideal_degree_rows(G, 6)
print("  deg 6 exact dim_Q =", dim_exact(rows, 6), f"[{time.time()-t0:.0f}s]")
