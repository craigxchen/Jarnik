"""Coherence test (round 14): is every K-good system (agreeing with the genuine labels at all levels of G) coherent on S,
i.e. is there ONE conjugation pattern eps with the label in the eps-class {u (eps o L) mod n} at every q in S?
ILP: columns = distinct (signature, eps-mask) per modulus (eps-mask = set of eps-classes containing the label),
level-1 relaxation outside G; for every eps: some q in S picks a column whose mask misses eps.
INFEASIBLE => coherence on S (exact, given G).  Usage: python3 coh.py K "ps" "G" "S" """
import sys, os, math, ctypes, time, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
from rlib14 import *
from forced import forced_set
from deform import U1_orbit_reps
lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "sigc3.so")); lib.run.restype = ctypes.c_int64
def cols_at(ps, K, q, V, last, epss):
    r = len(ps); m, L = genuine_labels(ps, q); n = m // 4; nuv = nu(ps, q)
    feas = np.zeros((len(V), 4), np.uint8)
    for i, v in enumerate(V):
        nz = math.prod(p ** abs(x) for p, x in zip(ps, v))
        for s in range(4): feas[i, s] = 1 if q * q <= nz * (2 if s % 2 else 1) else 0
    units = [u for u in range(1, n) if math.gcd(u, n) == 1] if n > 1 else [0]
    allc = []; off = [0]
    for e in epss:
        cl = sorted({sum(((u * e[j] * L[j]) % n) * n ** (r - 1 - j) for j in range(r)) for u in units}) if n > 1 else [0]
        allc += cl; off.append(len(allc))
    cls = np.array(allc, np.int64); offa = np.array(off, np.int32)
    reps, _ = U1_orbit_reps(m, nuv[0]); nw = (len(V) * 4 + 63) // 64 + 1; cap = 1 << 16
    Varr = np.ascontiguousarray(np.array(V, np.int32).reshape(-1)); fe = np.ascontiguousarray(feas.reshape(-1))
    nuarr = np.array(nuv, np.int32); rep = np.array(reps, np.int32)
    while True:
        tab = np.zeros(cap * nw, np.uint64); flg = np.zeros(cap, np.uint8)
        P = lambda a: ctypes.c_void_p(a.ctypes.data)
        cnt = lib.run(r, m, len(V), P(Varr), P(fe), P(last), P(cls), len(epss), P(offa), P(nuarr), len(reps), P(rep),
                      ctypes.c_int64(cap), P(tab), P(flg))
        if cnt >= 0: break
        cap *= 4
    bymask = {}
    for pos in np.nonzero(flg)[0]:
        words = tab[pos * nw:(pos + 1) * nw].tolist(); mask = int(words[-1]); x = 0
        for w in reversed(words[:-1]): x = (x << 64) | int(w)
        bymask.setdefault(mask, set()).add(x)
    out = []
    for mask, sigs in bymask.items():
        lst = list(sigs)
        if len(lst) <= 12000:
            mins = []
            for s in sorted(lst, key=lambda y: bin(y).count("1")):
                if not any((k & s) == k for k in mins): mins.append(s)
            lst = mins
        out += [(s, mask) for s in lst]
    return out
if __name__ == "__main__":
    K = int(sys.argv[1]); ps = [int(x) for x in sys.argv[2].split(",")]
    G = [int(x) for x in sys.argv[3].split(",")] if sys.argv[3] else []
    S = [int(x) for x in sys.argv[4].split(",")]
    r = len(ps); epss = [e for e in itertools.product([1, -1], repeat=r) if e[0] == 1]
    F = forced_set(ps, K, 60000); V = VK_mod_sign(r, K); pis = gaussian_primes(ps)
    last = np.array([max(j for j in range(r) if v[j]) for v in V], np.int32)
    Fr = [q for q in F if q not in G]; cols = []
    t0 = time.time()
    for q in Fr:
        cs = cols_at(ps, K, q, V, last, epss); cols += [(q, s, mk) for s, mk in cs]
    print(f"columns: {len(cols)} ({time.time()-t0:.0f}s)", flush=True)
    B = []; cons = []
    for v in V:
        z = Zv(pis, v); w = sum(abs(x) * math.log(p) for x, p in zip(v, ps))
        for s in range(4):
            B.append(w / 2 + (math.log(2) / 2 if s % 2 else 0.0)); cons.append(sum(vq(Ds(z, s), q) * math.log(q) for q in G))
    P_ = len(V) * 4; NE = len(epss); qi = {q: k for k, q in enumerate(Fr)}
    A = lil_matrix((P_ + len(Fr) + NE, len(cols)))
    for c, (q, s, mk) in enumerate(cols):
        x = s; b = 0; lq = math.log(q)
        while x:
            if x & 1: A[b, c] = lq
            x >>= 1; b += 1
        A[P_ + qi[q], c] = 1
        if q in S:
            for k in range(NE):
                if not (mk >> k) & 1: A[P_ + len(Fr) + k, c] = 1
    lb = np.concatenate([np.full(P_, -np.inf), np.ones(len(Fr)), np.ones(NE)])
    ub = np.concatenate([np.array(B) - np.array(cons) + 1e-9, np.ones(len(Fr)), np.full(NE, np.inf)])
    res = milp(c=np.zeros(len(cols)), constraints=LinearConstraint(csr_matrix(A), lb, ub), integrality=np.ones(len(cols)),
               bounds=Bounds(0, 1), options=dict(time_limit=7200))
    print(f"K={K} G={G} S={S}: status {res.status}: {'INFEASIBLE -> every good system is coherent on S' if res.status == 2 else res.message} ({time.time()-t0:.0f}s)", flush=True)
    if res.status == 0:
        ch = {cols[c][0]: cols[c][2] for c in range(len(cols)) if res.x[c] > 0.5}
        print("  eps-masks chosen on S:", {q: ch[q] for q in S})
