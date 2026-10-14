# Independent spanning certificate: "monomials of degree <= D span all functions on A" (mod p => over Q),
# via isotypic decomposition under a cyclic coordinate permutation sigma (order q) times k -> -k.
import numpy as np, flint, itertools, sys, time
def is_prime(n):
    if n < 2: return False
    i = 2
    while i*i <= n:
        if n % i == 0: return False
        i += 1
    return True
def find_prime(q, below=2**20):
    p = below - (below % (2*q)) + 1
    while not is_prime(p) or p >= below: p -= 2*q
    return p
def prim_root_of_unity(q, p):
    for g in range(2, p):
        z = pow(g, (p-1)//q, p)
        if all(pow(z, q//r, p) != 1 for r in range(2, q+1) if q % r == 0 and is_prime(r)) and z != 1:
            return z
from math import comb
import itertools as _it
def stratified(d, D, eps, need, perm, q, rng):
    """distinct sigma-orbit representatives of exponent vectors in {0,1,2}^d, degree <= D, parity eps,
    spread over degrees: equal quota per degree, all orbits of a degree taken when there are few."""
    degs = [t for t in range(D + 1) if t % 2 == eps]
    cnt = {t: sum(comb(d, j) * comb(d - j, t - 2*j) for j in range(0, t//2 + 1) if t - 2*j <= d - j) for t in degs}
    def canon(e):
        ee = e.copy(); best = tuple(ee)
        for _ in range(q - 1):
            ee = ee[perm]; best = min(best, tuple(ee))
        return best
    quota = {}; remaining = need; open_ = list(degs)
    while open_:
        share = -(-remaining // len(open_))
        small = [t for t in open_ if cnt[t] // q + 1 <= share]
        if not small:
            for t in open_: quota[t] = share
            break
        for t in small:
            quota[t] = cnt[t]; remaining -= cnt[t] // q + 1; open_.remove(t)
    mons = []
    for t in degs:
        seen = set(); want = quota[t]
        if cnt[t] <= 20000:
            for twos in range(0, t//2 + 1):
                ones = t - 2*twos
                if ones + twos > d: continue
                for P2 in _it.combinations(range(d), twos):
                    rest = [x for x in range(d) if x not in P2]
                    for P1 in _it.combinations(rest, ones):
                        e = np.zeros(d, dtype=np.int64); e[list(P2)] = 2; e[list(P1)] = 1
                        k = canon(e)
                        if k not in seen: seen.add(k); mons.append(np.array(k))
        else:
            w = np.array([comb(d, j) * comb(d - j, t - 2*j) if t - 2*j <= d - j else 0 for j in range(t//2 + 1)], float)
            tries = 0
            while len(seen) < want and tries < 50 * want:
                tries += 1
                j = rng.choice(len(w), p=w / w.sum())
                pos = rng.permutation(d); e = np.zeros(d, dtype=np.int64)
                e[pos[:j]] = 2; e[pos[j:j + t - 2*j]] = 1
                k = canon(e)
                if k not in seen: seen.add(k); mons.append(np.array(k))
    return mons
def certify(A, perm, D, q, label, seed=1, extra=300, chis=None):
    """A: int8 array N x d, invariant under (sigma k)_i = k_{perm[i]} and negation."""
    t0 = time.time()
    N, d = A.shape
    p = find_prime(q); z = prim_root_of_unity(q, p)
    pw3 = 3**np.arange(d, dtype=np.int64)
    code = ((A.astype(np.int64)+1) * pw3).sum(1)
    order = np.argsort(code); A = A[order]; code = code[order]
    assert np.all(np.diff(code) > 0)
    def idx_of(X):
        c = ((X.astype(np.int64)+1) * pw3).sum(1); i = np.searchsorted(code, c)
        assert np.all(code[i] == c), "A not invariant"
        return i
    imgs = [np.arange(N)]; X = A
    for s in range(1, q):
        X = X[:, perm]; imgs.append(idx_of(X))
    assert np.all(idx_of(A[:, perm] if False else X[:, perm]) == np.arange(N)), "sigma^q != id"
    imgs = np.stack(imgs, 1)              # N x q : index of sigma^s k
    neg = idx_of(-A)
    allimg = np.concatenate([imgs, neg[imgs]], 1)   # images under (s,+),(s,-)
    canon = allimg.min(1)
    reps = np.unique(canon)
    R = imgs[reps]                        # reps x q
    stabp = (R == reps[:, None])          # sigma^s k = k
    stabm = (neg[R] == reps[:, None])     # sigma^s k = -k
    # monomial evaluation via bitmasks
    nz = ((A != 0).astype(np.int64) * (1 << np.arange(d))).sum(1)
    ng = ((A < 0).astype(np.int64) * (1 << np.arange(d))).sum(1)
    par = np.array([bin(i).count("1") & 1 for i in range(1 << d)], dtype=np.int8)
    rng = np.random.default_rng(seed)
    results = {}
    chis = chis if chis is not None else [(j, e) for e in (0, 1) for j in range(q)]
    for (j, eps) in chis:
        valid = np.array([all((j*s) % q == 0 for s in np.nonzero(stabp[i])[0]) and
                          all((j*s) % q == 0 and eps == 0 for s in np.nonzero(stabm[i])[0]) for i in range(len(reps))])
        rows = int(valid.sum())
        # monomials of degree <= D with degree parity eps, distinct sigma-orbits
        mons = stratified(d, D, eps, rows + extra, perm, q, rng)
        w = np.array([pow(z, (-(j*s)) % q, p) for s in range(q)], dtype=np.float64)
        Rv = R[valid]
        Bm = np.zeros((rows, len(mons)), dtype=np.int64)
        for c0 in range(0, len(mons), 64):
            batch = mons[c0:c0+64]
            supp = np.array([int(((e > 0) * (1 << np.arange(d))).sum()) for e in batch])
            odd = np.array([int(((e % 2 == 1) * (1 << np.arange(d))).sum()) for e in batch])
            vals = np.where((nz[:, None] & supp[None, :]) == supp[None, :],
                            1 - 2*par[ng[:, None] & odd[None, :]].astype(np.int64), 0)   # N x b
            G = vals[Rv]                       # rows x q x b
            Bm[:, c0:c0+len(batch)] = (np.einsum('rsb,s->rb', G.astype(np.float64), w) % p).astype(np.int64)
        M = flint.nmod_mat(rows, len(mons), p)
        for i_ in range(rows):
            rw = Bm[i_].tolist()
            for j_, v_ in enumerate(rw):
                if v_: M[i_, j_] = v_
        del Bm
        r = M.rank()
        results[(j, eps)] = (rows, r)
        print(f"{label} chi=(j={j},eps={eps}) rows={rows} cols={len(mons)} rank={r} {'FULL' if r == rows else 'DEFICIENT'} t={time.time()-t0:.0f}s", flush=True)
    tot = sum(r for r, _ in results.values())
    print(label, "p", p, "sum of rows", tot, "N", N, "all full:", all(a == b for a, b in results.values()), flush=True)
    return results
