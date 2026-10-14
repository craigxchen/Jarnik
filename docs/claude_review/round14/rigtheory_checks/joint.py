"""Exact joint-rigidity test in the level-1 relaxation (round 14).
Usage: python3 joint.py K "ps" q0 [q0 ...]      (q0 = 0: every forced modulus)
For each q in F_K^ex: all distinct level-1 signatures {(v,s)} of parity-respecting labels whose collisions are each
individually affordable (log q <= B(v,s)); only inclusion-minimal ones are kept (a subset signature never costs more).
ILP (scipy/HiGHS): one signature per modulus, sum over moduli of log q over signatures containing (v,s) <= B(v,s),
and at q0 a signature realised by a label outside the genuine class.  Infeasible => every K-good system (all levels)
on F_K^ex is in the genuine class at q0 (exact joint rigidity; the relaxation only enlarges the feasible set).
Feasible => a candidate (level-1) joint deformation, to be lifted (lift.py)."""
import sys, os, math, ctypes, time, pickle
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
from rlib14 import *
from forced import forced_set
from deform import U1_orbit_reps
EPS = os.environ.get("EPS", "1") == "1"
HERE = os.path.dirname(os.path.abspath(__file__))
lib = ctypes.CDLL(os.path.join(HERE, "sigc.so"))
lib.run.restype = ctypes.c_int64
def signatures(ps, K, q, V, last):
    r = len(ps); pis = gaussian_primes(ps)
    m, L = genuine_labels(ps, q); n = m // 4; nuv = nu(ps, q)
    feas = np.zeros((len(V), 4), np.uint8)
    for i, v in enumerate(V):
        nz = math.prod(p ** abs(x) for p, x in zip(ps, v))
        for s in range(4): feas[i, s] = 1 if q * q <= nz * (2 if s % 2 else 1) else 0
    units = [u for u in range(1, n) if math.gcd(u, n) == 1] if n > 1 else [0]
    import itertools
    epss = [e for e in itertools.product([1, -1], repeat=r) if e[0] == 1] if EPS else [tuple([1] * r)]
    cls = np.array(sorted({sum(((u * e[j] * L[j]) % n) * n ** (r - 1 - j) for j in range(r)) for u in units for e in epss}) if n > 1 else [0], np.int64)
    reps, _ = U1_orbit_reps(m, nuv[0])
    nw = (len(V) * 4 + 63) // 64
    cap = 1 << 16
    while True:
        tab = np.zeros(cap * nw, np.uint64); flg = np.zeros(cap, np.uint8)
        Varr = np.ascontiguousarray(np.array(V, np.int32).reshape(-1)); fe = np.ascontiguousarray(feas.reshape(-1))
        nuarr = np.array(nuv, np.int32); rep = np.array(reps, np.int32)
        P = lambda a: ctypes.c_void_p(a.ctypes.data)
        cnt = lib.run(r, m, len(V), P(Varr), P(fe), P(last), P(cls), len(cls), P(nuarr), len(reps), P(rep),
                      ctypes.c_int64(cap), P(tab), P(flg))
        if cnt >= 0: break
        cap *= 4
    sigs = []
    for pos in np.nonzero(flg)[0]:
        words = tab[pos * nw:(pos + 1) * nw]
        x = 0
        for w in reversed(words.tolist()): x = (x << 64) | int(w)
        sigs.append((x, int(flg[pos]) & 3))
    gsig = 0
    for i, v in enumerate(V):
        t = sum(a * b for a, b in zip(v, L)) % m
        if t % n == 0: gsig |= 1 << (i * 4 + t // n)
    return sigs, gsig, n
def minimal(sigset):
    if len(sigset) > 12000:   # skip the quadratic filter: non-minimal columns are dominated, harmless, only slower
        return list(sigset)
    out = []
    for s in sorted(sigset, key=lambda x: bin(x).count("1")):
        if not any((k & s) == k for k in out): out.append(s)
    return out
def build(ps, K, cache=None):
    F = forced_set(ps, K, 60000)
    V = VK_mod_sign(len(ps), K)
    last = np.array([max(j for j in range(len(ps)) if v[j]) for v in V], np.int32)
    data = {}
    QMAX = int(os.environ.get("QMAX", "10000000"))
    for q in F:
        t0 = time.time()
        if q > QMAX:
            data[q] = dict(allmin=[], nonclass=[], gsig=0, n=(q - chi4(q)) // 4, nsig=0); continue
        sigs, gsig, n = signatures(ps, K, q, V, last)
        allmin = minimal({s for s, f in sigs})
        nonclass = minimal({s for s, f in sigs if f & 2})
        data[q] = dict(allmin=allmin, nonclass=nonclass, gsig=gsig, n=n, nsig=len(sigs))
        print(f"  q={q} n={n}: distinct feasible signatures={len(sigs)}, minimal={len(allmin)}, minimal non-class={len(nonclass)} ({time.time()-t0:.1f}s)", flush=True)
    return F, V, data
def solve(ps, K, F, V, data, q0, extra_fixed=None):
    B = []
    for v in V:
        w = sum(abs(x) * math.log(p) for x, p in zip(v, ps))
        for s in range(4): B.append(w / 2 + (math.log(2) / 2 if s % 2 else 0.0))
    cols = []  # (q, sig)
    for q in F:
        lst = data[q]['nonclass'] if q == q0 else data[q]['allmin']
        for s in lst: cols.append((q, s))
    nvar = len(cols); P = len(V) * 4
    A = lil_matrix((P + len(F), nvar))
    qi = {q: k for k, q in enumerate(F)}
    for c, (q, s) in enumerate(cols):
        lq = math.log(q); x = s; b = 0
        while x:
            if x & 1: A[b, c] = lq
            x >>= 1; b += 1
        A[P + qi[q], c] = 1
    lb = np.concatenate([np.full(P, -np.inf), np.ones(len(F))]); ub = np.concatenate([np.array(B) + 1e-9, np.ones(len(F))])
    res = milp(c=np.zeros(nvar), constraints=LinearConstraint(csr_matrix(A), lb, ub), integrality=np.ones(nvar),
               bounds=Bounds(0, 1), options=dict(time_limit=3600))
    sol = None
    if res.status == 0:
        sol = {cols[c][0]: cols[c][1] for c in range(nvar) if res.x[c] > 0.5}
    return res.status, res.message, sol, nvar
if __name__ == "__main__":
    K = int(sys.argv[1]); ps = [int(x) for x in sys.argv[2].split(",")]; q0s = [int(x) for x in sys.argv[3:]]
    F, V, data = build(ps, K)
    pickle.dump((F, V, data), open(f"joint_data{'_eps' if EPS else ''}{os.environ.get('QTAG','')}_K{K}_{'_'.join(map(str, ps))}.pkl", "wb"))
    if q0s == [0]: q0s = F
    for q0 in q0s:
        if not data[q0]['nonclass']:
            print(f"K={K} q0={q0}: no individually-feasible non-class label at all -> rigid (trivially)", flush=True); continue
        t0 = time.time()
        st, msg, sol, nvar = solve(ps, K, F, V, data, q0)
        if st == 2 or 'nfeasible' in msg:
            print(f"K={K} q0={q0}: ILP INFEASIBLE ({nvar} vars, {time.time()-t0:.1f}s) -> joint-rigid at q0", flush=True)
        elif st == 0:
            dev = [q for q in F if sol[q] != data[q]['gsig']]
            print(f"K={K} q0={q0}: FEASIBLE ({nvar} vars, {time.time()-t0:.1f}s); signature differs from genuine at {dev}", flush=True)
        else:
            print(f"K={K} q0={q0}: solver status {st}: {msg}", flush=True)
