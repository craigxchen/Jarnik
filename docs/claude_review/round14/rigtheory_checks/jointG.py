"""Joint deformation with a fixed genuine set G, and exact all-level witness construction (round 14).
Usage: python3 jointG.py K "ps" "G" q0 [q0 ...]       (G comma list, may be empty "")
ILP as in joint.py over q in F \\ G (level-1 minimal signatures); q in G carry genuine labels at ALL levels and
consume the genuine q-part of D_s(Z_v).  If feasible: (i) reconstruct an actual parity-respecting vector at every
q in F \\ G realising the chosen signature (outside the genuine class at q0), (ii) lift each to level 2 with NO level-2
collision of any v in V_K (random search; then no deeper collisions exist), (iii) verify every budget exactly in
integers: m_(v,s)^2 <= Norm(Z_v) * (2 if s odd).  A verified witness is a K-good system at all levels, on all of
F_K^ex, outside the genuine class at q0, so q0 is NOT in Core_J(K) (exact)."""
import sys, os, math, pickle, random, time
import numpy as np
from rlib14 import *
from forced import forced_set
from deform import U1_orbit_reps
from eps_class import eps_classes
EPS = os.environ.get("EPS", "1") == "1"
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
def ilp(ps, K, F, V, data, q0, G):
    pis = gaussian_primes(ps)
    B = []; cons = []
    for v in V:
        z = Zv(pis, v); w = sum(abs(x) * math.log(p) for x, p in zip(v, ps))
        for s in range(4):
            B.append(w / 2 + (math.log(2) / 2 if s % 2 else 0.0))
            D = Ds(z, s); cons.append(sum(vq(D, q) * math.log(q) for q in G))
    Fr = [q for q in F if q not in G]
    cols = []
    for q in Fr:
        for s in (data[q]['nonclass'] if q == q0 else data[q]['allmin']): cols.append((q, s))
    P = len(V) * 4; nvar = len(cols)
    A = lil_matrix((P + len(Fr), nvar)); qi = {q: k for k, q in enumerate(Fr)}
    for c, (q, s) in enumerate(cols):
        lq = math.log(q); x = s; b = 0
        while x:
            if x & 1: A[b, c] = lq
            x >>= 1; b += 1
        A[P + qi[q], c] = 1
    ub = np.concatenate([np.array(B) - np.array(cons) + 1e-9, np.ones(len(Fr))])
    lb = np.concatenate([np.full(P, -np.inf), np.ones(len(Fr))])
    res = milp(c=np.zeros(nvar), constraints=LinearConstraint(csr_matrix(A), lb, ub), integrality=np.ones(nvar),
               bounds=Bounds(0, 1), options=dict(time_limit=3600))
    if res.status != 0: return res.status, res.message, None
    return 0, res.message, {cols[c][0]: cols[c][1] for c in range(nvar) if res.x[c] > 0.5}
def find_vector(ps, q, V, target, need_nonclass):
    r = len(ps); m, L = genuine_labels(ps, q); n = m // 4; nuv = nu(ps, q)
    tv = np.zeros(len(V), bool); ts = np.zeros(len(V), np.int64)
    for i in range(len(V)):
        for s in range(4):
            if (target >> (i * 4 + s)) & 1: tv[i] = True; ts[i] = s
    units = [u for u in range(1, n) if math.gcd(u, n) == 1] if n > 1 else [0]
    import itertools
    epss = [e for e in itertools.product([1, -1], repeat=r) if e[0] == 1] if EPS else [tuple([1] * r)]
    cls = {tuple((u * e * x) % n for e, x in zip(ep, L)) for u in units for ep in epss} if n > 1 else None
    Varr = np.array(V, np.int64)
    reps, _ = U1_orbit_reps(m, nuv[0])
    ch = [np.arange(nuv[j], m, 2) for j in range(r)]
    for a in reps:
        rest = np.array(np.meshgrid(*ch[1:], indexing='ij')).reshape(r - 1, -1).T
        for k0 in range(0, len(rest), 200000):
            blk = rest[k0:k0 + 200000]
            grid = np.hstack([np.full((len(blk), 1), a), blk])
            t = (grid @ Varr.T) % m
            coll = (t % n == 0); s = (t // n) % 4
            ok = (coll == tv[None, :]).all(axis=1) & ((~coll) | (s == ts[None, :])).all(axis=1)
            for idx in np.nonzero(ok)[0]:
                ell = [int(x) for x in grid[idx]]
                inc = (tuple(x % n for x in ell) in cls) if cls is not None else True
                if need_nonclass and inc: continue
                return ell
    return None
def lift2(ps, q, V, ell, tries=200000, seed=1):
    r = len(ps); m = q - chi4(q); M2 = m * q; quarter = M2 // 4
    Varr = np.array(V, np.int64); rng = np.random.default_rng(seed)
    base = np.array(ell, np.int64)
    coll1 = [i for i, v in enumerate(V) if (sum(a * b for a, b in zip(v, ell)) % (m // 4)) == 0]
    if not coll1: return [int(x) for x in base], 0
    Vc = Varr[coll1]
    for k0 in range(0, tries, 5000):
        T = rng.integers(0, q, size=(5000, r))
        L2 = (base[None, :] + m * T) % M2
        t = (L2 @ Vc.T) % M2
        bad = (t % quarter == 0).any(axis=1)
        good = np.nonzero(~bad)[0]
        if len(good): return [int(x) for x in L2[good[0]]], 0
    return None, len(coll1)
def verify(ps, K, F, V, G, labels, lev2):
    """exact all-level check; labels[q] level-1 vectors for q not in G; lev2[q] level-2 lifts (collision-free)."""
    pis = gaussian_primes(ps); worst = None; viol = 0
    for v in V:
        z = Zv(pis, v); nz = z[0] ** 2 + z[1] ** 2
        for s in range(4):
            D = Ds(z, s); mm = 1
            for q in G: mm *= q ** vq(D, q)
            for q in F:
                if q in G: continue
                m = q - chi4(q); n = m // 4
                tt = sum(a * b for a, b in zip(v, labels[q])) % m
                if tt == s * n:
                    mm *= q
                    M2 = m * q; t2 = sum(a * b for a, b in zip(v, lev2[q])) % M2
                    assert t2 % (M2 // 4) != 0, "level-2 collision"
            if mm * mm > nz * (2 if s % 2 else 1): viol += 1
            ratio = math.log(mm) - 0.5 * math.log(nz * (2 if s % 2 else 1))
            worst = ratio if worst is None or ratio > worst else worst
    return viol, worst
if __name__ == "__main__":
    K = int(sys.argv[1]); ps = [int(x) for x in sys.argv[2].split(",")]
    G = [int(x) for x in sys.argv[3].split(",")] if sys.argv[3] else []
    q0s = [int(x) for x in sys.argv[4:]]
    F, V, data = pickle.load(open(f"joint_data{'_eps' if EPS else ''}{os.environ.get('QTAG','')}_K{K}_{'_'.join(map(str, ps))}.pkl", "rb"))
    for q0 in q0s:
        t0 = time.time()
        st, msg, sol = ilp(ps, K, F, V, data, q0, G)
        if st != 0:
            print(f"K={K} G={G} q0={q0}: ILP status {st} ({msg}) ({time.time()-t0:.0f}s)", flush=True); continue
        labels = {}; lev2 = {}; fail = None
        for q in F:
            if q in G: continue
            ell = find_vector(ps, q, V, sol[q], need_nonclass=(q == q0))
            if ell is None: fail = ("no vector", q); break
            labels[q] = ell
            l2, nb = lift2(ps, q, V, ell)
            if l2 is None: fail = ("no collision-free level-2 lift", q); break
            lev2[q] = l2
        if fail:
            print(f"K={K} G={G} q0={q0}: ILP feasible but witness construction failed: {fail}", flush=True); continue
        viol, worst = verify(ps, K, F, V, G, labels, lev2)
        inclass = bool(eps_classes(ps, q0, labels[q0]))
        nong = [q for q in F if q not in G and not eps_classes(ps, q, labels[q])]
        print(f"K={K} G={G} q0={q0}: WITNESS budget violations={viol} max log(m/budget)={worst:.4f} q0 in class={inclass} "
              f"moduli outside every eps-class={nong} ({time.time()-t0:.0f}s)", flush=True)
        pickle.dump(dict(K=K, ps=ps, G=G, q0=q0, labels=labels, lev2=lev2), open(f"witness{'_eps' if EPS else ''}_K{K}_q0_{q0}_G{'-'.join(map(str,G))}.pkl", "wb"))
