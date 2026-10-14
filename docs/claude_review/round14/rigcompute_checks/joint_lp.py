# joint_lp.py -- joint (everything-free) rigidity of a modulus q0, via a Lagrangian / LP-duality certificate.
#
# Joint problem J(q0, target): labels at EVERY forced q are free (free moduli: collision-free, zero demand).
# A good system exists with a target label at q0 iff one can pick, for every forced q, an achievable signature
# sigma_q (sigma_q0 achievable by a target label), with  sum_{q : (v,s) in sigma_q} log q <= B(v,s)  for all (v,s)
# (level-1 relaxation: deeper levels only add demand).  Signatures may be replaced by inclusion-minimal ones
# (demand is monotone), at q0 minimal among target-achievable ones.
# Certificate of infeasibility (weak duality): weights y_(v,s) >= 0 with
#      sum_q  log q * min_{sigma in Sigma_q} y(sigma)  >  sum_(v,s) y_(v,s) B(v,s),      y(sigma) = sum_{(v,s) in sigma} y_(v,s)
# (any feasible choice would give LHS <= sum_q log q * y(sigma_q) = sum y_(v,s) * demand <= RHS).
# The certificate is found by LP (HiGHS) and then VERIFIED EXACTLY: y rationalised, min over the complete
# exact signature list, logs enclosed in mpmath intervals.
# target: 'none' (label outside every column-sign class eps*L), 'notL' (outside class(L)).
import sys, pickle, math, time, json
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix, csr_matrix
from mpmath import iv, mp
from sigtab import minimal
from rl import cap

iv.prec = 120

def lp_cert(tabs, P, V, q0, target, Fset=None):
    Fset = Fset or sorted(tabs)
    fams = {}
    for q in Fset:
        T = tabs[q]
        if q == q0:
            fams[q] = minimal(T["sigs"], (lambda d: d["none"]) if target == "none" else (lambda d: d["notL"]))
        else:
            fams[q] = minimal(T["sigs"])
        if not fams[q]: return dict(status="no-target-label", q0=q0)
    pairs = sorted(set(e for q in Fset for d in fams[q] for e in d["sig"]))
    pid = {e: i for i, e in enumerate(pairs)}
    ny, nz = len(pairs), len(Fset)
    B = np.array([math.log(cap(P, V[i], s)) / 2 for (i, s) in pairs])
    rows = sum(len(fams[q]) for q in Fset)
    A = lil_matrix((rows, ny + nz)); k = 0
    for qi, q in enumerate(Fset):
        lq = math.log(q)
        for d in fams[q]:
            for e in d["sig"]: A[k, pid[e]] = -lq
            A[k, ny + qi] = 1.0; k += 1
    c = np.concatenate([B, -np.ones(nz)])           # minimise  sum y B - sum z
    bounds = [(0, 1)] * ny + [(None, None)] * nz
    res = linprog(c, A_ub=csr_matrix(A), b_ub=np.zeros(rows), bounds=bounds, method="highs")
    if res.status != 0: return dict(status="lp-fail", msg=res.message, q0=q0)
    gap = -res.fun
    out = dict(q0=q0, target=target, lp_gap=gap, status="lp-feasible-relaxation")
    if gap <= 1e-9: return out
    # exact verification with rationalised y
    y = [Fraction(max(0.0, v)).limit_denominator(10 ** 6) for v in res.x[:ny]]
    lhs = iv.mpf(0)
    for q in Fset:
        mn = min(sum((y[pid[e]] for e in d["sig"]), Fraction(0)) for d in fams[q])
        lhs += iv.log(q) * iv.mpf(mn.numerator) / mn.denominator
    rhs = iv.mpf(0)
    for (i, s), yy in zip(pairs, y):
        if yy: rhs += iv.log(cap(P, V[i], s)) / 2 * iv.mpf(yy.numerator) / yy.denominator
    ok = lhs.a > rhs.b
    out.update(status="CERTIFIED-INFEASIBLE" if ok else "lp-gap-not-verified", lhs=float(lhs.a), rhs=float(rhs.b),
               nsupport=sum(1 for v in y if v))
    if ok:
        out["cert"] = {str(pairs[i]): str(y[i]) for i in range(ny) if y[i]}
    return out

if __name__ == '__main__':
    D = pickle.load(open(sys.argv[1], "rb")); target = sys.argv[2]
    q0s = [int(x) for x in sys.argv[3].split(",")]
    from rl import monos_norm
    V = monos_norm(len(D["P"]), D["K"])
    for q0 in q0s:
        t0 = time.time(); o = lp_cert(D["tabs"], D["P"], V, q0, target)
        print({k: v for k, v in o.items() if k != "cert"}, f"[{time.time()-t0:.1f}s]", flush=True)
