"""LP relaxation of the joint2 ILP (round 14): LP infeasible => ILP infeasible => exact joint rigidity at q0.
Usage: python3 joint2_lp.py pickle K "ps" "G" q0 [q0 ...]"""
import sys, math, pickle, time
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix, csr_matrix
from rlib14 import *
from joint2 import coeffs
fn = sys.argv[1]; K = int(sys.argv[2]); ps = [int(x) for x in sys.argv[3].split(",")]
G = [int(x) for x in sys.argv[4].split(",")] if sys.argv[4] else []
F, V, data = pickle.load(open(fn, "rb"))
pis = gaussian_primes(ps); B = []; cons = []
for v in V:
    z = Zv(pis, v); w = sum(abs(x) * math.log(p) for x, p in zip(v, ps))
    for s in range(4):
        B.append(w / 2 + (math.log(2) / 2 if s % 2 else 0.0)); cons.append(sum(vq(Ds(z, s), q) * math.log(q) for q in G))
for q0 in [int(x) for x in sys.argv[5:]]:
    t0 = time.time(); Fr = [q for q in F if q not in G]; cols = []
    for q in Fr:
        for s in (data[q]['nonclass'] if q == q0 else data[q]['allmin']): cols.append((q, s))
    P = len(V) * 4; A = lil_matrix((P, len(cols))); Aeq = lil_matrix((len(Fr), len(cols))); qi = {q: k for k, q in enumerate(Fr)}
    for c, (q, s) in enumerate(cols):
        for row, k in coeffs(s, data[q]['lev'], len(V)).items(): A[row, c] = k * math.log(q)
        Aeq[qi[q], c] = 1
    res = linprog(np.zeros(len(cols)), A_ub=csr_matrix(A), b_ub=np.array(B) - np.array(cons) + 1e-9, A_eq=csr_matrix(Aeq),
                  b_eq=np.ones(len(Fr)), bounds=(0, 1), method="highs")
    print(f"K={K} G={G} q0={q0}: LP status {res.status} ({'INFEASIBLE -> joint-rigid' if res.status == 2 else res.message}) ({time.time()-t0:.0f}s, {len(cols)} vars)", flush=True)
