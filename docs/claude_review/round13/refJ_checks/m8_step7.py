"""All rational planes spanned by two low-width directions (the 358 primitive v with h_A(v) <= 7): do their pullbacks
(degree <= 6) fill the degree-6 gap J(A)_6 / (J_w)_6 (dimension 21)?  Planes deduplicated by their reduced row echelon form."""
import itertools, math, pickle, time
import numpy as np
import sympy as sp
from m8lib import *
A, cols, H, JA, W = pickle.load(open('m8_step1.pkl','rb'))
gw = pickle.load(open('m8_gw.pkl','rb'))
rows_w6 = ideal_degree_rows(gw, 6); base = dim_exact(rows_w6, 6)
dirs = [v for v, t, h in W]
seen = set(); t0 = time.time(); cache = {}
contrib = []; ncontrib = 0; nplanes = 0
for p, q in itertools.combinations(dirs, 2):
    M = sp.Matrix([p, q]).rref()[0]
    key = tuple(M)
    if key in seen: continue
    seen.add(key); nplanes += 1
    g, rho, ny = proj_ideal_gens(A, (p, q), 6)
    if 6 in g:
        ncontrib += 1
        contrib.append(((p, q), rho, [pullback(f, (p, q), cache) for f in g[6]]))
print(f"{nplanes} distinct planes from {len(dirs)} low-width directions; {ncontrib} have degree-6 generators of J(Y) [{time.time()-t0:.0f}s]", flush=True)
allpb = [f for _, _, l in contrib for f in l]
print("dim_Q (J_w + all these plane pullbacks)_6 =", dim_exact(rows_w6 + allpb, 6), " vs J_w:", base, " J(A)_6:", JA[6], f"[{time.time()-t0:.0f}s]")
rh = {}
for _, rho, _ in contrib: rh[rho] = rh.get(rho, 0) + 1
print("rho(Y) histogram (computed up to degree 6; None = rho > 6) of contributing planes:", rh)
pickle.dump(contrib, open('m8_step7.pkl', 'wb'))
