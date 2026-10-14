import numpy as np, itertools, time, sys
from math import comb
from inc_rank import Echelon, find_prime
from cyc_rank import prim_root_of_unity
def certify_inc(A, perm, D, q, label, chis, seed=1, batch=400, maxenum=4000):
    t0 = time.time(); N, d = A.shape
    p = find_prime(q); z = prim_root_of_unity(q, p)
    pw3 = 3**np.arange(d, dtype=np.int64)
    code = ((A.astype(np.int64)+1) * pw3).sum(1); order = np.argsort(code); A = A[order]; code = code[order]
    def idx_of(X):
        c = ((X.astype(np.int64)+1) * pw3).sum(1); i = np.searchsorted(code, c); assert np.all(code[i] == c); return i
    imgs = [np.arange(N)]; X = A
    for s in range(1, q): X = X[:, perm]; imgs.append(idx_of(X))
    assert np.all(idx_of(X[:, perm]) == np.arange(N))
    imgs = np.stack(imgs, 1); neg = idx_of(-A)
    canon = np.concatenate([imgs, neg[imgs]], 1).min(1); reps = np.unique(canon)
    R = imgs[reps]; stabp = (R == reps[:, None]); stabm = (neg[R] == reps[:, None])
    nz = ((A != 0).astype(np.int64) * (1 << np.arange(d))).sum(1)
    ng = ((A < 0).astype(np.int64) * (1 << np.arange(d))).sum(1)
    par = np.array([bin(i).count("1") & 1 for i in range(1 << d)], dtype=np.int8)
    rng = np.random.default_rng(seed); allfull = True; tot = 0
    def canon_e(e):
        ee = e.copy(); best = tuple(ee)
        for _ in range(q - 1): ee = ee[perm]; best = min(best, tuple(ee))
        return best
    for (j, eps) in chis:
        valid = np.array([all((j*s) % q == 0 for s in np.nonzero(stabp[i])[0]) and
                          all((j*s) % q == 0 and eps == 0 for s in np.nonzero(stabm[i])[0]) for i in range(len(reps))])
        rows = int(valid.sum()); Rv = R[valid]
        w = np.array([pow(z, (-(j*s)) % q, p) for s in range(q)], dtype=np.float64)
        E = Echelon(rows, p); log = []
        def cols(batch_e):
            supp = np.array([int(((e > 0) * (1 << np.arange(d))).sum()) for e in batch_e])
            odd = np.array([int(((e % 2 == 1) * (1 << np.arange(d))).sum()) for e in batch_e])
            out = np.zeros((rows, len(batch_e)))
            for s in range(q):
                ii = Rv[:, s]
                v = np.where((nz[ii][:, None] & supp[None, :]) == supp[None, :],
                             1 - 2*par[ng[ii][:, None] & odd[None, :]].astype(np.int64), 0)
                out += v * w[s]
            return out % p
        for t in [t for t in range(D + 1) if t % 2 == eps]:
            if len(E.piv) == rows: break
            cnt = sum(comb(d, a) * comb(d - a, t - 2*a) for a in range(t//2 + 1) if t - 2*a <= d - a)
            seen = set(); gained = 0
            if cnt <= maxenum:
                allm = []
                for a in range(t//2 + 1):
                    if t - 2*a > d - a: continue
                    for P2 in itertools.combinations(range(d), a):
                        rest = [x for x in range(d) if x not in P2]
                        for P1 in itertools.combinations(rest, t - 2*a):
                            e = np.zeros(d, dtype=np.int64); e[list(P2)] = 2; e[list(P1)] = 1
                            k = canon_e(e)
                            if k not in seen: seen.add(k); allm.append(np.array(k))
                for c0 in range(0, len(allm), batch):
                    gained += E.add(cols(allm[c0:c0+batch]))
            else:
                wts = np.array([comb(d, a) * comb(d - a, t - 2*a) if t - 2*a <= d - a else 0 for a in range(t//2 + 1)], float)
                tries = 0
                while True:
                    bm = []
                    while len(bm) < batch and tries < 40 * batch:
                        tries += 1
                        a = rng.choice(len(wts), p=wts / wts.sum()); pos = rng.permutation(d)
                        e = np.zeros(d, dtype=np.int64); e[pos[:a]] = 2; e[pos[a:a + t - 2*a]] = 1
                        k = canon_e(e)
                        if k not in seen: seen.add(k); bm.append(np.array(k))
                    if not bm: break
                    tries = 0
                    g = E.add(cols(bm)); gained += g
                    if g < 3 or len(E.piv) == rows: break
            log.append((t, gained))
        r = len(E.piv); tot += r; allfull &= (r == rows)
        print(f"{label} chi=(j={j},eps={eps}) rows={rows} rank={r} {'FULL' if r == rows else 'DEFICIENT'} per-degree gains={log} t={time.time()-t0:.0f}s", flush=True)
    print(label, "p", p, "sum rank", tot, "N", N, "all full", allfull, flush=True)
