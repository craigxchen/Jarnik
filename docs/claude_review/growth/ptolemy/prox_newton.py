r"""NUMERICAL exploration (not a proof): balanced degree-one rational curves in Mbar_{0,m+1}
via a proximal Gauss-Newton method in the X-chart (see balanced_newton.py for the chart).
Unknowns: coefficients of X_2..X_m (degree D=2^(m-2)-1; X_1=0; lead of X_2 fixed to 1) and
the roots rho_T of all blocks with 3<=|T|<=m-1 (two of them fixed for the affine gauge).
Each outer step solves  min ||F(z)||^2 + lam ||z_roots - target||^2  by Gauss-Newton and then
moves the target to the current roots while decreasing lam; this keeps the roots near a random
well-separated configuration and avoids the high-dimensional degenerate (merged-root) components.
Usage: python3 prox_newton.py m trials [seed]
"""
import sys, itertools
import numpy as np
from balanced_newton import setup, X, dX, pair_roots, verify

def run(m, rng, iters=400, verbose=False):
    rows, big, D, conds = setup(m)
    nb = len(big); E = len(conds); ncoef = (m - 1) * D + (m - 2)   # leads of X_3..X_m free
    c = lambda n: (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
    target = c(nb) * 3
    fixed = [0, 1]
    coef = {2: np.concatenate([[1.0], c(D)])}
    for i in rows[2:]:
        coef[i] = c(D + 1)
    rho = {T: target[n] for n, T in enumerate(big)}
    free = [T for n, T in enumerate(big) if n not in fixed]
    def pack():
        parts = [coef[2][1:]] + [coef[i] for i in rows[2:]]
        return np.concatenate(parts + [np.array([rho[T] for T in free])])
    def unpack(z):
        coef[2] = np.concatenate([[1.0], z[:D]])
        pos = D
        for i in rows[2:]:
            coef[i] = z[pos:pos + D + 1]; pos += D + 1
        for n, T in enumerate(free):
            rho[T] = z[pos + n]
    col = {}
    pos = 0
    col[2] = (0, D, 1)          # (start, length, offset into power vector)
    pos = D
    for i in rows[2:]:
        col[i] = (pos, D + 1, 0); pos += D + 1
    rootpos = {T: pos + n for n, T in enumerate(free)}
    nz = pos + len(free)
    def F():
        return np.array([X(coef, a, rho[T]) - X(coef, b, rho[T]) for T, a, b in conds])
    def J():
        Jm = np.zeros((E, nz), dtype=complex)
        for e, (T, a, b) in enumerate(conds):
            x = rho[T]
            pw = x ** np.arange(D, -1, -1)
            for i, s in ((a, 1), (b, -1)):
                if i == 1:
                    continue
                st, ln, off = col[i]
                Jm[e, st:st + ln] += s * pw[off:]
            if T in rootpos:
                Jm[e, rootpos[T]] = dX(coef, a, x) - dX(coef, b, x)
        return Jm
    z = pack()
    lam = 1.0
    tgt = np.array([rho[T] for T in free])
    for it in range(iters):
        unpack(z)
        f = F(); Jm = J()
        # augmented least squares: [J; sqrt(lam) S] dz = [-f; -sqrt(lam)(roots - tgt)]
        S = np.zeros((len(free), nz), dtype=complex)
        for n, T in enumerate(free):
            S[n, rootpos[T]] = 1.0
        rr = np.array([rho[T] for T in free]) - tgt
        A = np.vstack([Jm, np.sqrt(lam) * S])
        bvec = np.concatenate([-f, -np.sqrt(lam) * rr])
        dz = np.linalg.lstsq(A, bvec, rcond=None)[0]
        t = 1.0
        base = np.linalg.norm(f) ** 2 + lam * np.linalg.norm(rr) ** 2
        while t > 1e-8:
            zn = z + t * dz; unpack(zn)
            rrn = np.array([rho[T] for T in free]) - tgt
            if np.linalg.norm(F()) ** 2 + lam * np.linalg.norm(rrn) ** 2 < base:
                break
            t /= 2
        z = zn
        unpack(z)
        if it % 10 == 9:
            tgt = np.array([rho[T] for T in free]); lam *= 0.3
        if np.linalg.norm(F()) < 1e-13 and lam < 1e-8:
            break
    unpack(z)
    nr = np.linalg.norm(F())
    # Newton polish without proximal term
    for it in range(50):
        f = F()
        if np.linalg.norm(f) < 1e-14:
            break
        dz = np.linalg.lstsq(J(), -f, rcond=None)[0]
        z = z + dz; unpack(z)
    nr = np.linalg.norm(F())
    pr = pair_roots(m, rows, big, coef, rho)
    allrho = dict(rho); allrho.update(pr)
    vals = list(allrho.values())
    sep = min(abs(a - b) for a, b in itertools.combinations(vals, 2))
    leads = [0.0] + [coef[i][0] for i in rows[1:]]
    lsep = min(abs(x - y) for x, y in itertools.combinations(leads, 2))
    return nr, sep, lsep, coef, allrho, rows

if __name__ == "__main__":
    m = int(sys.argv[1]); trials = int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    rng = np.random.default_rng(seed)
    good = 0
    for tr in range(trials):
        nr, sep, lsep, coef, allrho, rows = run(m, rng)
        ok = nr < 1e-9 and sep > 1e-5 and lsep > 1e-5
        w = verify(m, rows, coef, allrho) if ok else float('nan')
        good += ok
        print(f"trial {tr}: residual {nr:.1e}, min sep {sep:.2e}, lead sep {lsep:.2e}, identity err {w:.1e} -> {'NONDEGENERATE' if ok else 'fail'}", flush=True)
    print("nondegenerate:", good, "/", trials)
