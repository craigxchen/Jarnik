"""Exact DFS (integer arithmetic) for the joint relaxation: one minimal signature per forced modulus (target at q0: a
minimal non-eps-class signature), every budget prod_q q^(a_q(v,s)) ^2 <= Norm(Z_v) * 2^[s odd].
G: moduli fixed to genuine all-level demand (as fixed products). Returns a feasible assignment or proves none."""
import sys, pickle, math, time
from glib import *
def canon(v):
    for x in v:
        if x: return x > 0
def run(ps, K, F, data, q0, G=(), fixed=None, limit=10**9):
    pis = [two_sq(p) for p in ps]; V = [v for v in VK(len(ps), K) if canon(v)]
    cap2 = []; P = []
    for v in V:
        z = Zv(pis, v); nz = z[0]**2 + z[1]**2
        for s in range(4):
            cap2.append(nz*(2 if s % 2 else 1)); g = 1
            for q in G: g *= q**vq(Ds(z, s), q)
            P.append(g)
    assert all(P[i]**2 <= cap2[i] for i in range(len(P)))
    mods = [q for q in F if q not in G]
    opts = {q: [tuple((row, q**a) for row, a in sig) for sig in (data[q]['nonclass'] if q == q0 else data[q]['allmin'])] for q in mods}
    nodes = [0]
    def ok(o):
        return all((P[row]*f)**2 <= cap2[row] for row, f in o)
    def rec(left):
        nodes[0] += 1
        if nodes[0] > limit: raise RuntimeError("limit")
        if not left: return {}
        best = None
        for q in left:
            c = [o for o in opts[q] if ok(o)]
            if not c: return None
            if best is None or len(c) < len(best[1]): best = (q, c)
        q, c = best; rest = [x for x in left if x != q]
        for o in c:
            for row, f in o: P[row] *= f
            r = rec(rest)
            for row, f in o: P[row] //= f
            if r is not None: r[q] = o; return r
        return None
    t0 = time.time(); r = rec(mods)
    return r, nodes[0], time.time()-t0
if __name__ == "__main__":
    ps = [int(x) for x in sys.argv[1].split(",")]; K = int(sys.argv[2]); fn = sys.argv[3]; q0 = int(sys.argv[4])
    G = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 and sys.argv[5] else []
    F, data = pickle.load(open(fn, "rb"))
    r, nodes, dt = run(ps, K, F, data, q0, G)
    print(f"K={K} q0={q0} G={G}: {'FEASIBLE' if r is not None else 'INFEASIBLE (exact)'} nodes={nodes} time={dt:.1f}s")
