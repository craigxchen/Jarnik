"""Find a small infeasible subset of moduli with HiGHS (integer-scaled conservative relaxation), then certify it by the
exact integer DFS of jdfs.py.  Dropping moduli only lowers demand, so an infeasible subset certifies the full problem."""
import sys, math, pickle, time, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
from glib import *
import jdfs
SC = 10**6
def canon(v):
    for x in v:
        if x: return x > 0
def ilp(ps, K, data, q0, mods, G=()):
    pis = [two_sq(p) for p in ps]; V = [v for v in VK(len(ps), K) if canon(v)]
    cap = []
    for v in V:
        z = Zv(pis, v); nz = z[0]**2 + z[1]**2
        for s in range(4):
            g = 1
            for q in G: g *= q**vq(Ds(z, s), q)
            B = 0.5*math.log(nz*(2 if s % 2 else 1)) - math.log(g)
            cap.append(math.ceil(B*SC) + 1)          # upper bound on budget (relaxation)
    cols = []
    for q in mods:
        for sig in (data[q]['nonclass'] if q == q0 else data[q]['allmin']): cols.append((q, sig))
    A = lil_matrix((len(cap) + len(mods), len(cols))); qi = {q: i for i, q in enumerate(mods)}
    for c, (q, sig) in enumerate(cols):
        for row, a in sig: A[row, c] = a*math.floor(math.log(q)*SC)   # lower bound on demand
        A[len(cap) + qi[q], c] = 1
    lb = np.concatenate([np.full(len(cap), -np.inf), np.ones(len(mods))]); ub = np.concatenate([np.array(cap, float), np.ones(len(mods))])
    res = milp(c=np.zeros(len(cols)), constraints=LinearConstraint(csr_matrix(A), lb, ub), integrality=np.ones(len(cols)), bounds=Bounds(0, 1),
               options=dict(time_limit=1800))
    return res.status   # 0 feasible, 2 infeasible
if __name__ == "__main__":
    ps = [int(x) for x in sys.argv[1].split(",")]; K = int(sys.argv[2]); fn = sys.argv[3]; q0 = int(sys.argv[4])
    G = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 and sys.argv[5] else []
    F, data = pickle.load(open(fn, "rb"))
    mods = [q for q in F if q not in G]
    t0 = time.time(); st = ilp(ps, K, data, q0, mods, G); print("full ILP status", st, f"{time.time()-t0:.0f}s", flush=True)
    if st != 2: sys.exit()
    for q in sorted(mods, reverse=True):
        if q == q0: continue
        trial = [x for x in mods if x != q]
        if ilp(ps, K, data, q0, trial, G) == 2: mods = trial
    print("infeasible subset (HiGHS):", mods, flush=True)
    r, nodes, dt = jdfs.run(ps, K, mods, data, q0, G)
    print(f"exact DFS on subset: {'FEASIBLE (subset)' if r is not None else 'INFEASIBLE (exact certificate)'} nodes={nodes} {dt:.1f}s", flush=True)
