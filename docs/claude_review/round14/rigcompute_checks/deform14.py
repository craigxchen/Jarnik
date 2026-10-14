# deform14.py -- one-modulus deformations of the genuine system on the exact forced set F_K^ex (round 14).
# System: genuine labels (all levels) at every q in F \ {q0}; collision-free labels at every free modulus (zero
# demand); an arbitrary parity-respecting ell at q0, level-1 relaxation (each level-1 collision at q0 costs log q0).
# For each q0: exhaustive DFS (./csp) for level-1 survivors OUTSIDE the genuine class c = u L mod n.
#   RIGID   : DFS complete, no non-genuine survivor -> no good system that is genuine off q0 is non-genuine at q0
#             (relaxation: exact certificate at all levels).
#   FLEX    : explicit non-genuine witness, re-verified here independently (exact integer budgets), lifted to level
#             2 with no level-2 collision -> an explicit good system (all levels) non-genuine at q0.
#   VACUOUS : the genuine class contains every parity-respecting c mod n (rigidity carries no information).
# Usage: python3 deform14.py THRJSON [mode=wit|count] [q0list|all] [maxnodes] [qmax]
import sys, os, json, time, math, random
import numpy as np
from rl import *

class Setup:
    def __init__(self, P, K, F):
        self.P = P; self.K = K; self.F = sorted(F); self.r = r = len(P)
        self.V = monos_norm(r, K); self.Va = np.array(self.V, dtype=np.int64)
        self.pis = [gauss_prime(p) for p in P]
        self.Z = [Zv(self.pis, v) for v in self.V]
        self.cap = [[cap(P, v, s) for s in range(4)] for v in self.V]
        # D_s(Z_v) as python ints and as int64 (bounded check)
        self.D = [[Ds(z, s) for s in range(4)] for z in self.Z]
        mx = max(abs(d) for row in self.D for d in row)
        assert mx < 2 ** 62
        Da = np.array(self.D, dtype=np.int64)
        # genuine F-smooth part of D_s(Z_v): Gfac[i][s] = prod_{q in F} q^{v_q(D)}; and per-q valuations
        self.val = {}          # q -> dict {(i,s): a}
        G = [[1] * 4 for _ in self.V]
        for q in self.F:
            idx = np.argwhere(Da % q == 0)
            dq = {}
            for i, s in idx:
                a = vq(self.D[i][s], q); dq[(int(i), int(s))] = a; G[i][s] *= q ** a
            self.val[q] = dq
        self.G = G
        for i in range(len(self.V)):
            for s in range(4):
                assert G[i][s] ** 2 <= self.cap[i][s], "genuine system not good?!"
    def other(self, q0):
        """G_other[i][s] = genuine demand product over F minus q0."""
        dq = self.val.get(q0, {})
        return [[self.G[i][s] // (q0 ** dq.get((i, s), 0)) for s in range(4)] for i in range(len(self.V))]

def modulus_data(S, q0):
    m = mq(q0); n = m // 4; T = Torus(q0)
    L = [T.L(pi) for pi in S.pis]; iota = T.iota
    # self-check: genuine level-1 collisions from logs == q0 | D_s(Z_v)
    dots = (S.Va @ np.array(L, dtype=np.int64)) % m
    for i in range(len(S.V)):
        for s in range(4):
            assert (dots[i] == (s * iota) % m) == (S.D[i][s] % q0 == 0)
    return m, n, L, iota

def verify_witness(S, q0, ell, Gother, m, n, iota, rng, tries=2000):
    """independent exact check of a level-1 witness + level-2 lift.  Returns (ok, collisions, lift)."""
    r = S.r; assert all(ell[j] % 2 == p for j, p in enumerate(nu(S.P, q0)))
    coll = []
    for i, v in enumerate(S.V):
        d = sum(a * b for a, b in zip(v, ell)) % m
        if d % n == 0:
            s = [t for t in range(4) if (t * iota) % m == d]; assert len(s) == 1; s = s[0]
            if (q0 * Gother[i][s]) ** 2 > S.cap[i][s]: return False, None, None
            coll.append((i, s))
    for _ in range(tries):
        t = [rng.randrange(q0) for _ in range(r)]
        ell2 = [ell[j] + m * t[j] for j in range(r)]
        if all(sum(a * b for a, b in zip(S.V[i], ell2)) % (n * q0) != 0 for i, s in coll):
            return True, coll, ell2
    return False, coll, None

def sign_classes(Lc, n, r):
    """union over column sign vectors eps (mod global sign) of {u * eps*L mod n}."""
    out = {}
    for bits in range(2 ** (r - 1)):
        eps = [1] + [(-1 if (bits >> (j - 1)) & 1 else 1) for j in range(1, r)]
        out[tuple(eps)] = genuine_class([(e * x) % n for e, x in zip(eps, Lc)], n)
    return out

def run_q0(S, q0, mode="wit", maxnodes=0, maxwit=3, rng=None, count_class=True, Shi=None):
    rng = rng or random.Random(q0)
    m, n, L, iota = modulus_data(S, q0)
    Go = S.other(q0)
    lq_ok = [[(q0 * Go[i][s]) ** 2 <= S.cap[i][s] for s in range(4)] for i in range(len(S.V))]
    forb = [[(s * iota) % m for s in range(4) if not lq_ok[i][s]] for i in range(len(S.V))]
    par = nu(S.P, q0)
    col0 = orbit_reps_U1(m, par[0])
    Lc = [x % n for x in L]
    cls = genuine_class(Lc, n)
    # vacuous? (class == all parity-respecting c mod n)
    npar = (n // 2) ** S.r if n % 2 == 0 else n ** S.r
    vac = (len(cls) == npar)
    res, wits = run_csp(m, n, S.r, 0 if mode == "count" else 1, maxwit, maxnodes, par, col0, S.V, forb, sorted(cls))
    out = dict(q0=q0, n=n, m=m, L=L, iota=iota, vacuous=vac, nallowed=sum(sum(x) for x in lq_ok), **res)
    ver = []
    for w in wits:
        Gv = Shi.other(q0) if Shi is not None else Go
        ok, coll, ell2 = verify_witness(Shi or S, q0, list(w), Gv, m, n, iota, rng)
        c = tuple(x % n for x in w)
        assert c not in cls
        ver.append(dict(ell=list(w), ok=ok, ncoll=len(coll) if coll else None,
                        coll=[(S.V[i], s) for i, s in coll] if coll else None, lift=ell2))
    out["witnesses"] = ver
    if res["status"] != "complete":
        out["verdict"] = "FLEX" if any(x["ok"] for x in ver) else "UNDECIDED"
    else:
        if res["nongen"] == 0 and not wits: out["verdict"] = "VACUOUS" if vac else "RIGID"
        else: out["verdict"] = "FLEX" if any(x["ok"] for x in ver) else "FLEX-L1-ONLY"
    return out

if __name__ == '__main__':
    thr = json.load(open(sys.argv[1])); P, K, F = thr["P"], thr["K"], thr["forced"]
    mode = sys.argv[2] if len(sys.argv) > 2 else "wit"
    qsel = sys.argv[3] if len(sys.argv) > 3 else "all"
    maxnodes = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    qmax = int(sys.argv[5]) if len(sys.argv) > 5 else 10 ** 9
    t0 = time.time(); S = Setup(P, K, F)
    ab = thr.get('aborted') or []
    # aborted (status unknown) moduli: zero demand for RIGID certificates (relaxation), genuine labels for witnesses
    Shi = Setup(P, K, sorted(F + ab)) if ab else None
    print(f"# P={P} K={K} |V|={len(S.V)} |F|={len(F)} q*={max(F)} setup {time.time()-t0:.1f}s", flush=True)
    q0s = F if qsel == "all" else [int(x) for x in qsel.split(",")]
    outs = []
    for q0 in q0s:
        if q0 > qmax: break
        t1 = time.time(); o = run_q0(S, q0, mode, maxnodes, maxwit=(20 if ab else 3), Shi=Shi)
        w = o["witnesses"][0] if o["witnesses"] else None
        print(f"q0={q0} n={o['n']} {o['verdict']} nodes={o['nodes']} survivors={o['leaves'] if mode=='count' else '-'} "
              f"nongen={o['nongen'] if mode=='count' else len(o['witnesses'])} allowed(v,s)={o['nallowed']} "
              f"wit={w['ell'] if w else None} ncoll={w['ncoll'] if w else None} [{time.time()-t1:.1f}s]", flush=True)
        outs.append(o)
    fn = sys.argv[1].replace("thr_", f"def_{mode}_").replace(".json", ".json")
    json.dump(outs, open(fn, "w"), default=str)
