"""Numerical exploration (NOT a proof): polynomial (function-field) solutions of the norm-level
Ptolemy/matching system with linear blocks, for m nonanchor rows (k=m+1 points).

Model: rows X_1..X_m in C[s] (or R[s]), X_m = 0 after translation, Y_i = 1.  Every block
n_T = s - r_T is linear (T subset [m], 2<=|T|<=m-1; the full block is a common factor and is
dropped).  Requirement: for every pair i<j,
      X_i - X_j = const * prod_{T contains i,j} (s - r_T),
with all roots r_T pairwise distinct.  Equivalently X_i(r_T) = X_j(r_T) for i,j in T, |T|>=3,
plus a nondegeneracy check that the remaining (doubleton) roots are new and distinct.

Parametrize X_i = (lam_i + mu_i s) Q_i(s), Q_i = prod_{T containing i and m, 3<=|T|<=m-1} (s-r_T).
Unknowns: some roots + (lam,mu); square Newton system after fixing random values for the rest.
Usage: python3 poly_family_newton.py m trials [real]
"""
import sys, itertools
import numpy as np

m = int(sys.argv[1]); trials = int(sys.argv[2]); REAL = len(sys.argv) > 3
rng = np.random.default_rng(12345 + m)

rows = list(range(1, m + 1))
mblocks = [frozenset(T) for r in range(3, m) for T in itertools.combinations(rows, r) if m in T]
fblocks = [frozenset(T) for r in range(3, m) for T in itertools.combinations(rows[:-1], r)]
allblocks = mblocks + fblocks
conds = []
for T in fblocks:
    Ts = sorted(T)
    for a, b in zip(Ts, Ts[1:]):
        conds.append((T, a, b))
E = len(conds); nrow = m - 1
F = E + 1 - 2 * nrow
print(f"m={m}: m-blocks {len(mblocks)}, free blocks {len(fblocks)}, conditions {E}, free roots needed {F}, total roots {len(allblocks)}")
assert F <= len(allblocks) - 2

def Xval(i, s, roots, c):
    lam, mu = c[2 * (i - 1)], c[2 * (i - 1) + 1]
    q = lam + mu * s
    for T in mblocks:
        if i in T:
            q = q * (s - roots[T])
    return q

def residual(zfree, freeidx, fixed, v0):
    roots = dict(fixed)
    for n, T in enumerate(freeidx):
        roots[T] = zfree[n]
    c = zfree[len(freeidx):]
    res = [Xval(a, roots[T], roots, c) - Xval(b, roots[T], roots, c) for T, a, b in conds]
    res.append(np.dot(c, v0) - 1.0)
    return np.array(res), roots, c

def newton(z, freeidx, fixed, v0, iters=200):
    for it in range(iters):
        r, _, _ = residual(z, freeidx, fixed, v0)
        nr = np.linalg.norm(r)
        if nr < 1e-13:
            return z, nr
        J = np.zeros((len(r), len(z)), dtype=z.dtype)
        h = 1e-7
        for k in range(len(z)):
            dz = np.zeros_like(z); dz[k] = h
            J[:, k] = (residual(z + dz, freeidx, fixed, v0)[0] - r) / h
        try:
            step = np.linalg.lstsq(J, -r, rcond=None)[0]
        except np.linalg.LinAlgError:
            return z, nr
        # damping
        t = 1.0
        while t > 1e-4:
            zn = z + t * step
            if np.linalg.norm(residual(zn, freeidx, fixed, v0)[0]) < nr:
                break
            t /= 2
        z = zn
    return z, np.linalg.norm(residual(z, freeidx, fixed, v0)[0])

def nondegenerate(roots, c):
    """Check exact Boolean root structure of every difference; return min root separation."""
    allr = dict(roots)
    lam = c[0::2]; mu = c[1::2]
    if np.min(np.abs(mu)) < 1e-6:
        return None
    # doubletons with row m: root of lam_i + mu_i s
    for i in rows[:-1]:
        allr[frozenset((i, m))] = -lam[i - 1] / mu[i - 1]
    # doubletons i<j<m: divide X_i - X_j by known roots
    S = np.linspace(-1, 1, 3)
    for i, j in itertools.combinations(rows[:-1], 2):
        if abs(mu[i - 1] - mu[j - 1]) < 1e-6:
            return None
        known = [allr[T] for T in allblocks if i in T and j in T]
        # X_i - X_j is degree 2^(m-2)-1; leading coeff mu_i - mu_j; find remaining root via values
        deg = 2 ** (m - 2) - 1
        assert len(known) == deg - 1
        # evaluate ratio at a generic point to get the last root
        s0 = 0.37 + 0.11j
        val = Xval(i, s0, roots, c) - Xval(j, s0, roots, c)
        pk = (mu[i - 1] - mu[j - 1]) * np.prod([s0 - r for r in known])
        r_last = s0 - val / pk
        # verify with a second point
        s1 = -0.53 + 0.29j
        val1 = Xval(i, s1, roots, c) - Xval(j, s1, roots, c)
        pred1 = (mu[i - 1] - mu[j - 1]) * np.prod([s1 - r for r in known]) * (s1 - r_last)
        if abs(val1 - pred1) > 1e-6 * (1 + abs(val1)):
            return None
        allr[frozenset((i, j))] = r_last
    vals = list(allr.values())
    sep = min(abs(a - b) for a, b in itertools.combinations(vals, 2))
    return sep, allr

found = 0
for tr in range(trials):
    dtype = float if REAL else complex
    def rnd(n):
        return rng.normal(size=n) if REAL else rng.normal(size=n) + 1j * rng.normal(size=n)
    perm = list(allblocks)
    rng.shuffle(perm)
    freeidx = perm[:F]
    fixed = {T: v for T, v in zip(perm[F:], rnd(len(perm) - F))}
    v0 = rnd(2 * nrow)
    z0 = np.concatenate([rnd(F), rnd(2 * nrow)]).astype(dtype)
    z, nr = newton(z0, freeidx, fixed, v0)
    if nr < 1e-11:
        r, roots, c = residual(z, freeidx, fixed, v0)
        nd = nondegenerate(roots, c)
        if nd is not None and nd[0] > 1e-4:
            found += 1
            print(f"trial {tr}: converged, residual {nr:.2e}, min separation of all {len(nd[1])} block roots {nd[0]:.3e}")
            if found >= 3:
                break
        else:
            print(f"trial {tr}: converged to degenerate solution (sep {None if nd is None else nd[0]})")
    else:
        print(f"trial {tr}: no convergence ({nr:.2e})")
print("nondegenerate solutions found:", found)
