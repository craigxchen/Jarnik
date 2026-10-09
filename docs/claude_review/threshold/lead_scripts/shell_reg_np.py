# General (non-symmetric) measures: max order of vanishing at the point 1 of Laurent polys whose
# exponents lie in Q_{s,t} = {k in Z^M: sum k = s, ||k||_1 <= t}  (equivalently the shell, since
# per-shell decoupling; we test the full ball slice which contains every shell <= t).
# order = smallest m such that polys of degree <= m in (k_1..k_{M-1}) interpolate the point set.
import numpy as np, itertools, sys
p = 2147483647
def pts(M, s, t, shell_only):
    out = []
    def rec(prefix, rem):
        j = len(prefix)
        if j == M - 1:
            last = s - sum(prefix)
            l1 = t - rem + abs(last)
            if (abs(last) <= rem) and (not shell_only or abs(last) == rem):
                out.append(tuple(prefix) + (last,))
            return
        for a in range(-rem, rem + 1):
            rec(prefix + [a], rem - abs(a))
    rec([], t)
    return out
def rank_mod(A):
    A = A.copy() % p; r = 0; rows, cols = A.shape
    for c in range(cols):
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: A[[r, i]] = A[[i, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        col = A[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            A[nzr] = (A[nzr] - (col[nzr, None] * A[r][None, :]) % p) % p
        r += 1
        if r == rows: break
    return r
def order(M, s, t, shell_only=True):
    S = pts(M, s, t, shell_only)
    n = len(S)
    X = np.array([k[:-1] for k in S], dtype=np.int64) % p
    m = 0
    cols = [np.ones(n, dtype=np.int64)]  # monomials as columns, built incrementally by degree
    frontier = {(): np.ones(n, dtype=np.int64)}
    allcols = [np.ones(n, dtype=np.int64)]
    while True:
        A = np.array(allcols, dtype=np.int64)  # (#mons) x n
        if rank_mod(A) == n: return m, n
        m += 1
        newf = {}
        for key, v in frontier.items():
            start = key[-1] if key else 0
            for j in range(start, M - 1):
                nk = key + (j,)
                newf[nk] = (v * X[:, j]) % p
        frontier = newf
        allcols.extend(newf.values())
if __name__ == '__main__':
    M = int(sys.argv[1]); tmax = int(sys.argv[2])
    for t in range(1, tmax + 1):
        for s in range(0, min(t, 3) + 1):
            if (t - s) % 2: continue
            m, n = order(M, s, t)
            print(f"M={M} t={t} s={s} |S|={n} order={m} ratio={m/t:.4f}", flush=True)
