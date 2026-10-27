# loc.py -- LOCAL (intrinsic) rigidity of a modulus q, valid for EVERY good system (all other moduli free).
# A level-1 collision (v,s) at q is possible in a good system only if q^2 <= cap(v,s) = Norm(Z_v)*(2 if s odd)
# (the collision alone costs log q).  'Admissible' labels at q: parity-respecting ell all of whose collisions with
# ||v||_1 <= K are of that kind.  Decided exactly by exhaustive DFS (./csp):
#   LOC  : every admissible label is in class(L)            LOCe : every admissible label is in some class(eps*L)
# K-independence: if p_min^(K+1) >= q^2, no monomial of norm > K can give an intrinsic constraint at q (every such v
# has Norm(Z_v) >= p_min^(K+1) >= q^2), so the verdict is the same for every K' >= K ('stable').
# Usage: python3 loc.py NAME P K qmax [maxnodes]
import sys, json, time, math
from rl import *
from sigtab import class_masks

def loc_q(P, V, caps, pis, q, maxnodes=0):
    m = mq(q); n = m // 4; T = Torus(q); L = [T.L(pi) for pi in pis]; iota = T.iota
    forb = [[(s * iota) % m for s in range(4) if q * q > caps[i][s]] for i in range(len(V))]
    par = nu(P, q); col0 = orbit_reps_U1(m, par[0])
    out = {}
    for tgt in ("eps", "L"):
        cls = set(class_masks(L, n, len(P)).keys()) if tgt == "eps" else genuine_class([x % n for x in L], n)
        res, w = run_csp(m, n, len(P), 1, 1, maxnodes, par, col0, V, forb, sorted(cls))
        if res["status"] != "complete" and not w: out[tgt] = None
        else: out[tgt] = (len(w) == 0, list(w[0]) if w else None)
        if tgt == "eps" and out[tgt] and not out[tgt][0]:
            out["L"] = (False, out[tgt][1]); break
    return n, out

if __name__ == '__main__':
    name = sys.argv[1]; P = [int(x) for x in sys.argv[2].split(",")]; K = int(sys.argv[3]); qmax = int(sys.argv[4])
    maxnodes = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    r = len(P); V = monos_norm(r, K); pis = [gauss_prime(p) for p in P]
    caps = [[cap(P, v, s) for s in range(4)] for v in V]
    rows = []
    for q in primerange(3, qmax + 1):
        if q in P: continue
        t0 = time.time(); n, o = loc_q(P, V, caps, pis, q, maxnodes)
        stable = min(P) ** (K + 1) >= q * q
        vac = n <= 2
        rows.append(dict(q=q, n=n, LOCe=o["eps"][0] if o["eps"] else None, LOC=o["L"][0] if o.get("L") else None,
                         wit=o["eps"][1] if o["eps"] else None, stable=stable))
        print(f"q={q} n={n} LOCe={rows[-1]['LOCe']} LOC={rows[-1]['LOC']} stable={stable} wit={rows[-1]['wit']} "
              f"[{time.time()-t0:.1f}s]", flush=True)
    json.dump(dict(P=P, K=K, rows=rows), open(f"loc_{name}_K{K}.json", "w"))
