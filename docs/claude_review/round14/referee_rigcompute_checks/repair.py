# repair.py -- try to turn an author's joint witness (jwmix pkl) into one that also satisfies the per-v angle
# compatibility (exact), by replacing the label at ONE forced modulus q0 (n_q0 >= 3), keeping non-coherence on
# F ∩ [3, Q].  New label at q0: level-1 label + random lift with no level-2 collision (integer-label semantics).
import sys, json, math, random, itertools, pickle
import numpy as np
from reflib import *
from rverify import R

def demands(P, K, labs, V, skip=None):
    dem = {}
    for q, top in labs.items():
        if q == skip: continue
        m = m_of(q); n = m // 4
        L, iota, _ = genuine_L(P, q); k = (iota // n) % 4
        for i, v in enumerate(V):
            d = sum(a * b for a, b in zip(v, top))
            if d % n: continue
            s = ((d % m) // n) * (1 if k == 1 else 3) % 4
            dem[(i, s)] = dem.get((i, s), 1) * q ** (1 + val(d // n, q))
    return dem

def main(fn, Ffile, maxq):
    W = R(open(fn, "rb")).load(); P, K, Q = W["P"], W["K"], W["Q"]; r = len(P)
    labs = {q: [int(x) for x in W["labels"][q]["lifts"][-1]] for q in W["labels"]}
    F = sorted(labs); V = monos(r, K)
    capi = [[capv(P, v, s) for s in range(4)] for v in V]; Nv = [norm_v(P, v) for v in V]
    full = (1 << 2 ** (r - 1)) - 1
    masks = {}
    for q in F:
        n = m_of(q) // 4; L, _, _ = genuine_L(P, q); masks[q] = class_mask([x % n for x in labs[q]], L, n, r)
    random.seed(7)
    for q0 in [q for q in F if m_of(q) // 4 >= 3 and q <= maxq]:
        base = demands(P, K, labs, V, skip=q0)
        m = m_of(q0); n = m // 4; L, iota, _ = genuine_L(P, q0); k = (iota // n) % 4
        accother = full
        for q in F:
            if q <= Q and q != q0: accother &= masks[q]
        nu = nu_par(P, q0)
        Vm = np.array(V, dtype=np.int64)
        cols = [np.arange(nu[j], m, 2, dtype=np.int64) for j in range(r)]
        allv = np.array(list(itertools.product(*cols)), dtype=np.int64)
        d = (allv @ Vm.T) % m; coll = (d % n) == 0; s = ((d // n) * (1 if k == 1 else 3)) % 4
        # per-unit precheck with base demand
        allowed = np.array([[(base.get((i, t), 1) * q0) ** 2 <= capi[i][t] for t in range(4)] for i in range(len(V))])
        ok = np.ones(len(allv), dtype=bool)
        for t in range(4): ok &= ~(coll & (s == t) & (~allowed[:, t][None, :])).any(axis=1)
        cand = np.nonzero(ok)[0]; random.shuffle(cand := list(cand))
        found = None
        for idx in cand[:20000]:
            ell = [int(x) for x in allv[idx]]
            mk = class_mask([x % n for x in ell], L, n, r)
            if q0 <= Q and (accother & mk) != 0: continue
            if q0 > Q and accother != 0: continue
            dem = dict(base)
            for i in np.nonzero(coll[idx])[0]:
                dem[(int(i), int(s[idx, i]))] = dem.get((int(i), int(s[idx, i])), 1) * q0
            if not all(angle_ok_exact([dem.get((i, t), 1) for t in range(4)], Nv[i]) for i in range(len(V))): continue
            # lift with no level-2 collision
            cv = [V[i] for i in np.nonzero(coll[idx])[0]]
            for _ in range(300):
                lift = [x + m * random.randrange(q0) for x in ell]
                if all(sum(a * b for a, b in zip(v, lift)) % (n * q0) != 0 for v in cv): break
            else: continue
            found = (q0, ell, lift); break
        print(fn.split("/")[-1], "q0", q0, "level-1 PU-ok candidates", len(cand), "repair:", found, flush=True)
        if found:
            labs2 = dict(labs); labs2[q0] = found[2]
            dem = demands(P, K, labs2, V)
            pu = all(dd * dd <= capi[i][t] for (i, t), dd in dem.items())
            ang = all(angle_ok_exact([dem.get((i, t), 1) for t in range(4)], Nv[i]) for i in range(len(V)))
            acc = full
            for q in F:
                if q <= Q:
                    n2 = m_of(q) // 4; L2, _, _ = genuine_L(P, q); acc &= class_mask([x % n2 for x in labs2[q]], L2, n2, r)
            print("   full re-check: per-unit", pu, "angle", ang, "non-coherent on F∩[3,Q]", acc == 0, flush=True)
            json.dump({str(q): v for q, v in labs2.items()}, open("repaired_" + fn.split("/")[-1].replace(".pkl", ".json"), "w"))
            return

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
