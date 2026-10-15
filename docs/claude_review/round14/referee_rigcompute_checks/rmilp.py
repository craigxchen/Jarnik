# rmilp.py -- referee's independent joint test (level-1 relaxation, Lemma R items 1,2,4,5), solved by HiGHS MILP.
# Moduli: forced q <= Qsub.  Per modulus: all individually admissible parity-respecting labels, grouped into
# (signature, column-sign mask) types.  Question: is there a level-1-good choice that is NOT coherent with any
# column-sign class on C = F ∩ [3, Q]?   INFEASIBLE (solver) = independent (floating-point) confirmation of Q_core>=Q.
import sys, json, math, itertools, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
from reflib import *

DEEP = {}
def types_at(P, K, q, V, capi):
    r = len(P); m1 = m_of(q); n = m1 // 4; a_ = DEEP.get(q, 1); m = m1 * q ** (a_ - 1)
    L, iota, _ = genuine_L(P, q); k = (iota // n) % 4
    nu = nu_par(P, q)
    cols = [np.arange(nu[j], m, 2, dtype=np.int64) for j in range(r)]
    U1 = [u for u in range(1, m, 4) if math.gcd(u, m) == 1]
    seen = set(); reps = []
    for x in cols[0]:
        x = int(x)
        if x in seen: continue
        orb = set((u * x) % m for u in U1); seen |= orb; reps.append(x)
    rest = np.array(list(itertools.product(*cols[1:])), dtype=np.int64)
    Vm = np.array(V, dtype=np.int64)
    adm = np.array([[q * q <= capi[i][s] for s in range(4)] for i in range(len(V))])
    T = {}
    cm = {}
    for kk, e in enumerate(eps_list(r)):
        for u in units(n):
            c = tuple((u * ei * x) % n for ei, x in zip(e, L)); cm[c] = cm.get(c, 0) | (1 << kk)
    if n == 1: cm = {tuple([0] * r): (1 << 2 ** (r - 1)) - 1}
    for x0 in reps:
        ell = np.concatenate([np.full((len(rest), 1), x0, dtype=np.int64), rest], axis=1)
        d = (ell @ Vm.T) % m
        coll = (d % n) == 0
        s = (((d % m1) // n) * (1 if k == 1 else 3)) % 4
        dep = np.where(coll, 1, 0)
        for b in range(1, a_):
            dep += (coll & ((d % (n * q ** b)) == 0)).astype(np.int64)
        ok = np.ones(len(ell), dtype=bool)
        for t in range(4):
            for b in range(1, a_ + 1):
                admb = np.array([[q ** (2 * b) <= capi[i][tt] for tt in range(4)] for i in range(len(V))])
                ok &= ~(coll & (s == t) & (dep == b) & (~admb[:, t][None, :])).any(axis=1)
        full = (1 << 2 ** (r - 1)) - 1
        for idx in np.nonzero(ok)[0]:
            sig = frozenset((int(i), int(s[idx, i]), int(dep[idx, i])) for i in np.nonzero(coll[idx])[0])
            ms = T.setdefault(sig, set())
            if full in ms: continue
            ms.add(cm.get(tuple(int(x) % n for x in ell[idx]), 0))
    out = [(sig, mk) for sig, ms in T.items() for mk in ms]
    return out

def main(name, K, Qsub, Q, Fsrc):
    thr = json.load(open(Fsrc)); P = thr["P"]; r = len(P)
    F = [q for q in sorted(thr["forced"]) if q <= Qsub]
    V = monos(r, K, half=True)
    capi = [[capv(P, v, s) for s in range(4)] for v in V]
    t0 = time.time()
    types = {q: types_at(P, K, q, V, capi) for q in F}
    print(name, K, "types per modulus:", {q: len(types[q]) for q in F}, f"[{time.time()-t0:.0f}s]", flush=True)
    idx = [];
    for q in F:
        for j, (sig, mask) in enumerate(types[q]): idx.append((q, sig, mask))
    nv = len(idx); nb = 2 ** (r - 1)
    rows = len(F) + len(V) * 4 + nb
    A = lil_matrix((rows, nv)); lo = np.zeros(rows); hi = np.zeros(rows)
    qpos = {q: i for i, q in enumerate(F)}
    for c, (q, sig, mask) in enumerate(idx):
        A[qpos[q], c] = 1
        for (i, s, b) in sig:
            A[len(F) + 4 * i + s, c] += b * math.log(q)
        if q <= Q:
            for b in range(nb):
                if not (mask >> b) & 1: A[len(F) + 4 * len(V) + b, c] = 1
    for i in range(len(F)): lo[i] = hi[i] = 1
    for i in range(len(V)):
        for s in range(4):
            rr = len(F) + 4 * i + s; lo[rr] = -np.inf; hi[rr] = 0.5 * math.log(capi[i][s])
    for b in range(nb):
        rr = len(F) + 4 * len(V) + b; lo[rr] = 1; hi[rr] = np.inf
    res = milp(c=np.zeros(nv), constraints=LinearConstraint(csr_matrix(A), lo, hi), integrality=np.ones(nv),
               bounds=Bounds(0, 1), options=dict(time_limit=4000, disp=False))
    print(name, K, "Qsub", Qsub, "Q", Q, "status", res.status, res.message, f"[{time.time()-t0:.0f}s]", flush=True)
    if res.x is not None:
        sel = [idx[c] for c in range(nv) if res.x[c] > 0.5]
        # exact re-check of level-1 budgets
        dem = {}
        for q, sig, mask in sel:
            for (i, s, b) in sig: dem[(i, s)] = dem.get((i, s), 1) * q ** b
        viol = [(V[i], s) for (i, s), d in dem.items() if d * d > capi[i][s]]
        acc = (1 << nb) - 1
        for q, sig, mask in sel:
            if q <= Q: acc &= mask
        print("  solution: exact level-1 violations", len(viol), "mask AND on C", acc, flush=True)
        json.dump([(q, sorted(sig), mask) for q, sig, mask in sel], open(f"rmilp_sol_{name}_K{K}_{Qsub}_{Q}.json", "w"))

if __name__ == "__main__":
    if len(sys.argv) > 6:
        for t in sys.argv[6].split(","):
            qq, aa = t.split(":"); DEEP[int(qq)] = int(aa)
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
