# One-modulus deformations of the genuine system (exact).
# System: genuine labels (all levels) at every forced q in F, collision-free above (Theorem 2).  Fix q0 in F and
# replace the label at q0 by an arbitrary parity-respecting ell.  Which ell keep the system good (norm <= K)?
# Enumerate ALL ell at q0 (as c in (Z/n)^r plus t in (Z/4)^r, ell = c + n t), exact per-unit level-1 test against
# the remaining budget; then for feasible ell outside the genuine class, an exact all-level lift search.
# Genuine class at q0 (level 1): c = u * L mod n with u a unit mod n.
import sys, math, itertools, time
import numpy as np
from lib import *
from forced_coh import load_set, Mod, lift_levels

def vq(x, q):
    x = abs(x); k = 0
    while x and x % q == 0: x //= q; k += 1
    return k

def run(K, r, tag, maxlift=40, only=None, extra=None, exclude=()):
    P, F = load_set(K, r, tag); V = monos(r, K); nv = len(V)
    B = np.array([[budget(P, v, s) for s in range(4)] for v in V])
    pis = [gauss_prime(p) for p in P]
    Zvs = [Zv(pis, v) for v in V]
    Gq = {q: np.array([[vq(Ds(z, s), q) * math.log(q) for s in range(4)] for z in Zvs]) for q in F}
    Gtot = sum(Gq[q] for q in F if q not in exclude)
    if extra is not None: Gtot = Gtot + extra
    assert (Gtot <= B + 1e-9).all()
    Va = np.array(V, dtype=np.int64)
    print(f"K={K} r={r} {tag} P={P} forced F={F}", flush=True)
    for q in F:
        if only and q not in only: continue
        t0 = time.time()
        md = Mod(q, P, K, V); m, n = md.m, md.n; lq = math.log(q)
        eps = 1 if md.iota == n else 3             # iota = eps * n
        slack = B - (Gtot - (Gq[q] if q not in exclude else 0))                 # remaining budget for the load at q (all levels)
        okS = slack >= lq - 1e-9                   # (v,s) can absorb one level-1 collision at q
        Lc = [x % n for x in md.L]
        units = [u for u in range(1, n + 1) if math.gcd(u, n) == 1] if n > 1 else [1]
        gen_class = set(tuple((u * x) % n for x in Lc) for u in units)
        ranges = [range(md.par[j], n, 2) if n % 2 == 0 else range(n) for j in range(r)]
        tot = feas = feas_gen = 0; nongen = []
        T4 = np.array(list(itertools.product(range(4), repeat=r)), dtype=np.int64)
        for chunk in itertools.islice(iter(lambda it=itertools.product(*ranges): list(itertools.islice(it, 200000)), []), 10**9):
            C = np.array(chunk, dtype=np.int64)
            D = C @ Va.T                                            # exact integers
            coll = (D % n == 0)
            d = (D // n) % 4
            # stage 1: each colliding v must have some absorbable unit
            anyok = okS.any(axis=1)
            pre = ~(coll & ~anyok[None, :]).any(axis=1)
            for ci in np.where(pre)[0]:
                c = C[ci]; cv = np.where(coll[ci])[0]
                # parity of ell_j = c_j + n t_j must be nu_j
                tmask = np.ones(len(T4), dtype=bool)
                if n % 2 == 1:
                    for j in range(r):
                        tmask &= ((c[j] + n * T4[:, j]) % 2 == md.par[j])
                Tt = T4[tmask]
                if len(cv) == 0:
                    tot += len(Tt); feas += len(Tt)
                    if tuple(int(x) for x in c) in gen_class: feas_gen += len(Tt)
                    else: nongen.extend([(c, t) for t in Tt[:2]])
                    continue
                svals = (eps * (d[ci, cv][None, :] + Tt @ Va[cv].T)) % 4     # unit index (iota = eps n)
                ok = okS[cv[None, :], svals].all(axis=1)
                tot += len(Tt); feas += int(ok.sum())
                if ok.any():
                    if tuple(int(x) for x in c) in gen_class: feas_gen += int(ok.sum())
                    elif len(nongen) < 5000: nongen.extend([(c, t) for t in Tt[ok][:2]])
            # count also the vectors rejected in stage 1 for 'tot'
            rej = np.where(~pre)[0]
            for ci in rej:
                tot += 4 ** r if n % 2 == 0 else 2 ** r
        # lift check for a sample of non-genuine feasible vectors
        lifted_ok = 0; checked = 0
        for (c, t) in nongen[:maxlift]:
            ell = [int(c[j] + n * t[j]) for j in range(r)]
            cl = md.collisions(ell)
            base = Gtot - (Gq[q] if q not in exclude else 0)
            nl = base.copy()
            for i, s in cl: nl[i, s] += lq
            if not (nl <= B + 1e-9).all(): continue
            checked += 1
            if lift_levels(md, ell, cl, nl, B, lq) is not None: lifted_ok += 1
        print(f"  q={q} n={n}: parity vectors at q {tot}; level-1 feasible {feas}; of which genuine class {feas_gen}; "
              f"non-genuine-class feasible {feas - feas_gen}; all-level lift found for {lifted_ok}/{checked} sampled non-genuine "
              f"[{time.time()-t0:.0f}s]", flush=True)

if __name__ == '__main__':
    K, r, tag = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    only = [int(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else None
    if len(sys.argv) > 5:      # pickle from forced_prefix: top moduli replaced by a saved non-genuine system
        import pickle
        pk = pickle.load(open(sys.argv[5], 'rb')); idx = int(sys.argv[6])
        sv = pk['saved'][idx]; top = sorted(sv['system'])
        print(f"top moduli {top} use saved non-genuine system #{idx}")
        run(K, r, tag, only=only, extra=sv['top_load'], exclude=set(top))
    else:
        run(K, r, tag, only=only)
