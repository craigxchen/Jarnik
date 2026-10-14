"""Joint-rigidity ILP with exact multi-level labels at chosen small primes (round 14).
Usage: [EPS=1] python3 joint2.py K "ps" "levels" "G" q0 [q0 ...]
  levels: e.g. "3:4,7:2" = labels at q=3 enumerated mod m_3*3^3 (levels 1..4 charged exactly), at 7 levels 1..2;
          every other forced q: level 1 only (higher levels relaxed, i.e. not charged).
  G     : moduli fixed to the genuine labels at ALL levels (consume the genuine q-part of D_s); "" for none.
  EPS=1 : the target at q0 is a level-1 label outside EVERY per-column-conjugation class {u (eps L)} (default);
          EPS=0 : outside the fixed class {u L} only.
The ILP is a relaxation of the true all-level problem (uncharged levels only lower demand), so INFEASIBLE means:
every K-good system (all levels) on F_K^ex that agrees with G is in the (eps-)genuine class at q0  (exact)."""
import sys, os, math, ctypes, time, pickle, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
from rlib14 import *
from forced import forced_set
from deform import U1_orbit_reps
EPS = os.environ.get("EPS", "1") == "1"
HERE = os.path.dirname(os.path.abspath(__file__))
lib = ctypes.CDLL(os.path.join(HERE, "sigc2.so")); lib.run.restype = ctypes.c_int64
def sigs_at(ps, K, q, V, last, lev):
    r = len(ps)
    Mtop, Ltop = genuine_labels(ps, q, a=lev)
    m1 = q - chi4(q); n = m1 // 4; nuv = nu(ps, q)
    L1 = [x % m1 for x in Ltop]
    feas = np.zeros((len(V), 4), np.uint8)
    for i, v in enumerate(V):
        nz = math.prod(p ** abs(x) for p, x in zip(ps, v))
        for s in range(4):
            d = 0
            while d < lev and q ** (2 * (d + 1)) <= nz * (2 if s % 2 else 1): d += 1
            feas[i, s] = d
    units = [u for u in range(1, n) if math.gcd(u, n) == 1] if n > 1 else [0]
    epss = [e for e in itertools.product([1, -1], repeat=r) if e[0] == 1] if EPS else [tuple([1] * r)]
    cls = np.array(sorted({sum(((u * e[j] * L1[j]) % n) * n ** (r - 1 - j) for j in range(r)) for u in units for e in epss})
                   if n > 1 else [0], np.int64)
    reps, _ = U1_orbit_reps(Mtop, nuv[0])
    nw = (len(V) * 4 * lev + 63) // 64; cap = 1 << 16
    while True:
        tab = np.zeros(cap * nw, np.uint64); flg = np.zeros(cap, np.uint8)
        Varr = np.ascontiguousarray(np.array(V, np.int32).reshape(-1)); fe = np.ascontiguousarray(feas.reshape(-1))
        nuarr = np.array(nuv, np.int32); rep = np.array(reps, np.int32)
        P = lambda a: ctypes.c_void_p(a.ctypes.data)
        cnt = lib.run(r, Mtop, len(V), P(Varr), P(fe), P(last), P(cls), len(cls), P(nuarr), len(reps), P(rep),
                      ctypes.c_int64(cap), P(tab), P(flg), q, lev)
        if cnt >= 0: break
        cap *= 4
    out = []
    for pos in np.nonzero(flg)[0]:
        x = 0
        for w in reversed(tab[pos * nw:(pos + 1) * nw].tolist()): x = (x << 64) | int(w)
        out.append((x, int(flg[pos]) & 3))
    return out
def minimal(sigset):
    out = []
    for s in sorted(sigset, key=lambda x: bin(x).count("1")):
        if not any((k & s) == k for k in out): out.append(s)
    return out
def build(ps, K, levels):
    F = forced_set(ps, K, 60000); V = VK_mod_sign(len(ps), K)
    last = np.array([max(j for j in range(len(ps)) if v[j]) for v in V], np.int32)
    data = {}
    for q in F:
        lev = levels.get(q, 1); t0 = time.time()
        sg = sigs_at(ps, K, q, V, last, lev)
        data[q] = dict(lev=lev, allmin=minimal({s for s, f in sg}), nonclass=minimal({s for s, f in sg if f & 2}), nsig=len(sg))
        print(f"  q={q} lev={lev}: signatures={len(sg)} minimal={len(data[q]['allmin'])} minimal-nonclass={len(data[q]['nonclass'])} ({time.time()-t0:.1f}s)", flush=True)
    return F, V, data
def coeffs(sig, lev, nV):
    """map signature bits (v,s,a) -> dict row(v*4+s) -> number of levels"""
    d = {}; x = sig; b = 0
    while x:
        if x & 1:
            row = b // lev; d[row] = d.get(row, 0) + 1
        x >>= 1; b += 1
    return d
def solve(ps, K, F, V, data, q0, G):
    pis = gaussian_primes(ps); B = []; cons = []
    for v in V:
        z = Zv(pis, v); w = sum(abs(x) * math.log(p) for x, p in zip(v, ps))
        for s in range(4):
            B.append(w / 2 + (math.log(2) / 2 if s % 2 else 0.0))
            cons.append(sum(vq(Ds(z, s), q) * math.log(q) for q in G))
    Fr = [q for q in F if q not in G]; cols = []
    for q in Fr:
        for s in (data[q]['nonclass'] if q == q0 else data[q]['allmin']): cols.append((q, s))
    P = len(V) * 4; A = lil_matrix((P + len(Fr), len(cols))); qi = {q: k for k, q in enumerate(Fr)}
    for c, (q, s) in enumerate(cols):
        for row, k in coeffs(s, data[q]['lev'], len(V)).items(): A[row, c] = k * math.log(q)
        A[P + qi[q], c] = 1
    ub = np.concatenate([np.array(B) - np.array(cons) + 1e-9, np.ones(len(Fr))]); lb = np.concatenate([np.full(P, -np.inf), np.ones(len(Fr))])
    res = milp(c=np.zeros(len(cols)), constraints=LinearConstraint(csr_matrix(A), lb, ub), integrality=np.ones(len(cols)),
               bounds=Bounds(0, 1), options=dict(time_limit=3600))
    sol = {cols[c][0]: cols[c][1] for c in range(len(cols)) if res.x is not None and res.x[c] > 0.5} if res.status == 0 else None
    return res.status, res.message, sol, len(cols)
if __name__ == "__main__":
    K = int(sys.argv[1]); ps = [int(x) for x in sys.argv[2].split(",")]
    levels = {int(a.split(":")[0]): int(a.split(":")[1]) for a in sys.argv[3].split(",") if a}
    G = [int(x) for x in sys.argv[4].split(",")] if sys.argv[4] else []
    q0s = [int(x) for x in sys.argv[5:]]
    fn = f"joint2_K{K}_{'_'.join(map(str, ps))}_lev{sys.argv[3].replace(':','-').replace(',','_')}_eps{int(EPS)}.pkl"
    if os.path.exists(fn): F, V, data = pickle.load(open(fn, "rb"))
    else:
        F, V, data = build(ps, K, levels); pickle.dump((F, V, data), open(fn, "wb"))
    if q0s == [0]: q0s = [q for q in F if q not in G]
    for q0 in q0s:
        if not data[q0]['nonclass']:
            print(f"K={K} lev={levels} G={G} eps={int(EPS)} q0={q0}: no individually feasible label outside the class -> rigid", flush=True); continue
        t0 = time.time(); st, msg, sol, nv = solve(ps, K, F, V, data, q0, G)
        tag = "INFEASIBLE -> joint-rigid at q0" if st == 2 else ("FEASIBLE (relaxation)" if st == 0 else f"status {st}: {msg}")
        print(f"K={K} lev={levels} G={G} eps={int(EPS)} q0={q0}: {tag} ({nv} vars, {time.time()-t0:.0f}s)", flush=True)
