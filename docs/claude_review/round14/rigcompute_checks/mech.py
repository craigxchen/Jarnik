# mech.py -- what makes a modulus rigid (round 14).  For each forced q0 of a set (one-modulus setting: genuine
# labels elsewhere), split the forbidden level-1 collisions (v,s) at q0 into
#   intrinsic : q0^2 > cap(v,s)                (no good system of any kind can collide there)
#   spent     : q0^2 <= cap but (q0 * G_other)^2 > cap   (forbidden only because genuine labels elsewhere used budget)
# and decide exactly (exhaustive DFS, ./csp):
#   LOC  : rigid with intrinsic constraints only  -> every admissible label at q0 lies in class(L)   [joint-valid]
#   LOCe : every admissible label at q0 lies in some column-sign class eps*L (target 'none' infeasible) [joint-valid,
#          up to per-modulus conjugation]
#   minimal norm k such that the constraints of monomials of norm <= k already force rigidity (one-modulus).
# Usage: python3 mech.py NAME K [qmax]
import sys, json, time, math
from rl import *
from deform14 import Setup, modulus_data
from sigtab import class_masks

def rigid_with(S, q0, m, n, L, iota, allowed_fn, target="L"):
    forb = [[(s * iota) % m for s in range(4) if not allowed_fn(i, s)] for i in range(len(S.V))]
    par = nu(S.P, q0); col0 = orbit_reps_U1(m, par[0])
    if target == "L":
        cls = genuine_class([x % n for x in L], n)
    else:   # union of all column-sign classes
        cls = set(class_masks(L, n, S.r).keys())
    res, wits = run_csp(m, n, S.r, 1, 1, 0, par, col0, S.V, forb, sorted(cls))
    assert res["status"] == "complete"
    return len(wits) == 0, (wits[0] if wits else None)

if __name__ == '__main__':
    name, K = sys.argv[1], int(sys.argv[2]); qmax = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
    thr = json.load(open(f"thr_{name}_K{K}.json")); P, F = thr["P"], thr["forced"]
    S = Setup(P, K, F); norms = [sum(abs(t) for t in v) for v in S.V]
    dw = {o["q0"]: o["verdict"] for o in json.load(open(f"def_wit_{name}_K{K}.json"))}
    rows = []
    for q0 in F:
        if q0 > qmax: break
        m, n, L, iota = modulus_data(S, q0); Go = S.other(q0)
        intr = lambda i, s: q0 * q0 <= S.cap[i][s]
        full = lambda i, s: (q0 * Go[i][s]) ** 2 <= S.cap[i][s]
        nint = sum(1 for i in range(len(S.V)) for s in range(4) if not intr(i, s))
        nsp = sum(1 for i in range(len(S.V)) for s in range(4) if intr(i, s) and not full(i, s))
        loc, _ = rigid_with(S, q0, m, n, L, iota, intr, "L")
        loce, _ = rigid_with(S, q0, m, n, L, iota, intr, "eps")
        kmin = None
        if dw.get(q0) in ("RIGID", "VACUOUS"):
            for k in range(1, K + 1):
                ok, _ = rigid_with(S, q0, m, n, L, iota, lambda i, s: full(i, s) or norms[i] > k, "L")
                if ok: kmin = k; break
        # genuine collisions at q0 and their norms
        gn = sorted(norms[i] for (i, s) in S.val.get(q0, {}))
        rows.append((q0, n, dw.get(q0), nint, nsp, loc, loce, kmin, len(gn), min(gn) if gn else None))
        print(f"q0={q0} n={n} onemod={dw.get(q0)} forbidden: intrinsic={nint} spent={nsp} LOC={loc} LOCe={loce} "
              f"k_min(one-mod)={kmin} genuine collisions={len(gn)} (min norm {min(gn) if gn else None})", flush=True)
    json.dump(rows, open(f"mech_{name}_K{K}.json", "w"))
