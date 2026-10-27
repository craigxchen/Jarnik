# indep_cert.py NAME K Q Qsub -- INDEPENDENT re-check of a GB coherence certificate (round 15).
# Own code only (numpy enumeration, own discrete logs, own class masks, own propagation); imports nothing from base/.
# Claim: every parity-respecting level-1 label choice on F ∩ [3,Qsub] with  prod_q q^[v collides at q with s] <= |D_s(Z_v)|
# for all v in V_K, s, is coherent (one column-sign class eps) on C = F ∩ [3,Q].  Since level-1 demands are <= the
# all-level demands and dropping moduli lowers demands, this certifies Q_core^GB(P,K) >= Q (round-14 Lemma R 2,5).
# Method: enumerate ALL parity-respecting ell in (Z/m)^r at each q (no symmetry reduction), keep the individually
# admissible ones (q <= |D_s| for each collision), then a sound fixpoint propagation:
#   mandatory(q, e) = log q if every remaining label at q collides with e;  drop a label at q if one of its collisions
#   e has  sum_(q' != q) mandatory(q', e) + log q  >  log|D_e| (+1e-9);  infeasible if some q has no label left;
#   coherence: the run is for NON-coherent systems, so if some eps lies in the class of every remaining label at every
#   q in C, no non-coherent system exists.  If the fixpoint does not decide, branch on the smallest domain (DFS).
import sys, json, math, itertools, time
import numpy as np
sys.setrecursionlimit(10000)

def gp(p):
    for a in range(1, int(math.isqrt(p)) + 1, 2):
        b2 = p - a * a; b = math.isqrt(b2)
        if b * b == b2 and b % 2 == 0: return (a, b)

def gm(x, y, q): return ((x[0] * y[0] - x[1] * y[1]) % q, (x[0] * y[1] + x[1] * y[0]) % q)

def torus(q):
    m = q - (1 if q % 4 == 1 else -1)
    elems = set()
    for a in range(q):
        for b in range(q):
            nn = (a * a + b * b) % q
            if nn: inv = pow(nn, -1, q); elems.add((((a * a - b * b) * inv) % q, (2 * a * b * inv) % q))
    assert len(elems) == m
    for g in sorted(elems):
        x = (1, 0); tab = {}
        for k in range(m):
            if x in tab: break
            tab[x] = k; x = gm(x, g, q)
        if len(tab) == m: return m, tab
    raise RuntimeError

def leg(a, q): return -1 if pow(a % q, (q - 1) // 2, q) == q - 1 else 1

def run(name, K, Q, Qsub, base_dir):
    thr = json.load(open(f"{base_dir}/thr_{name}_K{K}.json")); P, F = thr["P"], thr["forced"]
    r = len(P); pis = [gp(p) for p in P]
    V = [v for v in itertools.product(range(-K, K + 1), repeat=r) if 0 < sum(map(abs, v)) <= K and [t for t in v if t][-1] > 0]
    D = []
    for v in V:
        z = (1, 0)
        for pi, t in zip(pis, v):
            f = pi if t > 0 else (pi[0], -pi[1])
            for _ in range(abs(t)): z = (z[0] * f[0] - z[1] * f[1], z[0] * f[1] + z[1] * f[0])
        X, Y = z; D.append((Y, X - Y, X, X + Y))
    logD = {(i, s): math.log(abs(D[i][s])) for i in range(len(V)) for s in range(4)}
    Va = np.array(V, dtype=np.int64)
    qs = [q for q in F if q <= Qsub]; doms = {}; info = {}
    for q in qs:
        m, tab = torus(q); n = m // 4; iota = tab[(0, 1)]
        tor = lambda z: (lambda a, b: (((a * a - b * b) * pow((a * a + b * b) % q, -1, q)) % q, (2 * a * b * pow((a * a + b * b) % q, -1, q)) % q))(z[0] % q, z[1] % q)
        L = [tab[tor(pi)] for pi in pis]
        par = [1 if leg(p, q) == -1 else 0 for p in P]
        # sanity: genuine collisions <=> q | D_s
        for i, v in enumerate(V):
            d = sum(a * b for a, b in zip(v, L)) % m
            for s in range(4): assert (d == (s * iota) % m) == (D[i][s] % q == 0)
        grids = [np.arange(par[j], m, 2, dtype=np.int64) for j in range(r)]
        rest = np.array(np.meshgrid(*grids[1:], indexing="ij")).reshape(r - 1, -1).T
        forbid = [(i, (s * iota) % m) for i in range(len(V)) for s in range(4) if q > abs(D[i][s])]
        Es = []; Ds = []
        for x0 in grids[0]:
            E = np.hstack([np.full((len(rest), 1), x0, dtype=np.int64), rest])
            dots = (E @ Va.T) % m
            ok = np.ones(len(E), dtype=bool)
            for i, t in forbid: ok &= dots[:, i] != t
            Es.append(E[ok]); Ds.append(dots[ok])
        E = np.vstack(Es); dots = np.vstack(Ds)
        # class masks: bit b for eps vector (eps_0 = +1, bits for j >= 1)
        units = [u for u in range(1, n + 1) if math.gcd(u, n) == 1] if n > 1 else [1]
        cls = {}
        for b in range(2 ** (r - 1)):
            eps = [1] + [(-1 if (b >> (j - 1)) & 1 else 1) for j in range(1, r)]
            for u in units: cls.setdefault(tuple((u * e * x) % n for e, x in zip(eps, L)), set()).add(b)
        labs = []
        for k in range(len(E)):
            row = dots[k]; coll = frozenset((i, s) for i in np.nonzero(row % n == 0)[0] for s in range(4) if row[i] == (s * iota) % m)
            mask = 0
            for b in cls.get(tuple(int(x) % n for x in E[k]), ()): mask |= 1 << b
            labs.append((coll, mask))
        # merge identical (signature, mask)
        doms[q] = sorted(set(labs), key=lambda t: (len(t[0]), t[1])); info[q] = len(E)
    full = (1 << (2 ** (r - 1))) - 1; C = [q for q in qs if q <= Q]; lq = {q: math.log(q) for q in qs}
    nodes = [0]

    def propagate(dom):
        while True:
            mand = {}
            for q in qs:
                common = None
                for coll, mk in dom[q]: common = set(coll) if common is None else common & coll
                for e in common or (): mand[e] = mand.get(e, 0.0) + lq[q]
            for e, x in mand.items():
                if x > logD[e] + 1e-9: return False
            changed = False
            for q in qs:
                common = None
                for coll, mk in dom[q]: common = set(coll) if common is None else common & coll
                keep = []
                for coll, mk in dom[q]:
                    if all(mand.get(e, 0.0) - (lq[q] if e in common else 0.0) + lq[q] <= logD[e] + 1e-9 for e in coll): keep.append((coll, mk))
                if len(keep) < len(dom[q]): changed = True; dom[q] = keep
                if not keep: return False
            alive = 0
            for q in C:
                for coll, mk in dom[q]: alive |= (~mk) & full
            if alive != full: return False          # some eps common to every label at every q in C: coherent
            if not changed: return True

    def dfs(dom):
        nodes[0] += 1
        if not propagate(dom): return None
        open_ = [q for q in qs if len(dom[q]) > 1]
        if not open_:
            a = full
            for q in C: a &= dom[q][0][1]
            return None if a else dom
        qb = min(open_, key=lambda q: (q not in C, len(dom[q])))
        for opt in list(dom[qb]):
            d2 = dict(dom); d2[qb] = [opt]
            res = dfs(d2)
            if res: return res
        return None
    if VERB:
        d0 = {q: list(doms[q]) for q in qs}
        def masks(dd, q): return sorted(set(mk for _, mk in dd[q]))
        print("initial: q:#labels(signature,mask classes) and local status (LOCe = every admissible label in some eps-class)")
        for q in qs:
            loce = all(mk for _, mk in d0[q]); print(f"   q={q}: {len(d0[q])} distinct (sig,mask); LOCe={loce}; masks={[bin(x) for x in masks(d0, q)][:6]}")
        # propagate WITHOUT the coherence test to see forced reductions
        dd = {q: list(doms[q]) for q in qs}
        global_full = full
        ok = propagate_nc(dd) if 'propagate_nc' in globals() else None
    t0 = time.time(); res = dfs(dict(doms))
    print(f"{name} K={K} Q={Q} Qsub={Qsub}: admissible labels per q: " + " ".join(f"{q}:{info[q]}" for q in qs))
    print(f"RESULT {'NON-COHERENT SYSTEM FOUND (level 1)' if res else 'CERTIFIED: every GB level-1 system on F<=%d is coherent on F<=%d' % (Qsub, Q)} nodes={nodes[0]} [{time.time()-t0:.1f}s]")
    return res is None

VERB = False
if __name__ == '__main__':
    VERB = len(sys.argv) > 6
    run(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5] if len(sys.argv) > 5 else "base")
