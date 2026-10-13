"""Referee check of p3Dlower Theorem 2.1 (twin decomposition), exact modular arithmetic, independent code.
Random twin profiles (e = 1): a[:, f(x)] - a[:, s(x)] = eps_x e_x. Checks items 1-4 and the norm bounds, exhaustively
on small v, and that the decoupling reproduces any parity-compatible kappa."""
import random, itertools
import numpy as np
random.seed(11); ok = True
def profile(M, nJ):
    while True:
        eps = [random.choice((1, -1)) for _ in range(M)]
        J = np.random.randint(0, 2, size=(M, nJ)); S = np.random.randint(0, 2, size=(M, M))
        for x in range(M): S[x, x] = 0 if eps[x] == 1 else 1
        F = S.copy()
        for x in range(M): F[x, x] += eps[x]
        A = np.concatenate([J, F, S], axis=1)    # columns: J (0..nJ-1), f(x)=nJ+x, s(x)=nJ+M+x
        if len({tuple(r) for r in A}) == M: return A, eps
np.random.seed(5)
for trial in range(40):
    M = random.randint(2, 6); nJ = random.randint(1, 6); r = nJ + 2 * M
    A, eps = profile(M, nJ); m = random.choice([8, 12, 16, 24, 40, 44, 72])
    fcol = [nJ + x for x in range(M)]; scol = [nJ + M + x for x in range(M)]
    sig = A[:, scol]                           # sigma_{yx} = a_{y, s(x)}
    for _ in range(30):
        ell = np.random.randint(0, m, size=r)
        ellJ = ell[:nJ]; U = (ell[fcol] + ell[scol]) % m; Lam = (np.array(eps) * ell[fcol]) % m
        Fy = (A[:, :nJ] @ ellJ + sig @ U) % m
        kappa = (A @ ell) % m
        ok &= np.array_equal(kappa, (Fy + Lam) % m)                                   # item 2
        # inverse (item 1)
        lf = (np.array(eps) * Lam) % m; ls = (U - lf) % m
        ok &= np.array_equal(lf, ell[fcol] % m) and np.array_equal(ls, ell[scol] % m)
        # (2.2) for random v
        for _ in range(10):
            v = np.random.randint(-3, 4, size=r)
            cp = np.array(eps) * (v[fcol] - v[scol])
            lhs = int(v @ ell) % m; rhs = int(v[:nJ] @ ellJ + v[scol] @ U + cp @ Lam) % m
            ok &= (lhs == rhs)
        # item 4: decoupling. random parity vector nu, random parity-respecting (ellJ,U), kappa with right colouring
        nu = np.random.randint(0, 2, size=r)
        ellJ2 = (nu[:nJ] + 2 * np.random.randint(0, m // 2, size=nJ)) % m
        U2 = ((nu[fcol] + nu[scol]) % 2 + 2 * np.random.randint(0, m // 2, size=M)) % m
        F2 = (A[:, :nJ] @ ellJ2 + sig @ U2) % m
        colour = (A @ nu) % 2                                                   # sigma_q(y) = sum_j a_yj nu_j
        ok &= np.array_equal(colour, (F2 + nu[fcol]) % 2)                         # item 4 parity identity
        kap = (colour + 2 * np.random.randint(0, m // 2, size=M)) % m            # any kappa with the colouring
        Lam2 = (kap - F2) % m
        l = np.zeros(r, dtype=int); l[:nJ] = ellJ2
        l[fcol] = (np.array(eps) * Lam2) % m; l[scol] = (U2 - l[fcol]) % m
        ok &= np.array_equal(l % 2, nu) and np.array_equal((A @ l) % m, kap)
    # norm bounds and unimodularity (exhaustive small box, small r)
    if r <= 7:
        for v in itertools.product(range(-2, 3), repeat=r):
            v = np.array(v); cp = np.array(eps) * (v[fcol] - v[scol])
            vt = np.concatenate([v[:nJ], v[scol], cp]); n1, n2 = np.abs(v).sum(), np.abs(vt).sum()
            ok &= (n2 <= 2 * n1) and (n1 <= 2 * n2)
    Tm = np.zeros((r, r), dtype=int)                                            # matrix of v -> vt
    for i in range(nJ): Tm[i, i] = 1
    for x in range(M):
        Tm[nJ + x, scol[x]] = 1; Tm[nJ + M + x, fcol[x]] = eps[x]; Tm[nJ + M + x, scol[x]] = -eps[x]
    ok &= abs(round(np.linalg.det(Tm))) == 1
print("Thm 2.1 items 1-4, (2.2), norm bounds, unimodularity on 40 random twin profiles:", ok)
# support (sparsity) is NOT preserved: v = e_{s(x)} has 1-sparse v but 2-sparse vt
print("example: v = e_s(x) -> vt has support 2 (U_x and Lambda_x): 'same sparsity' holds only up to factor 2")
print("ALL OK" if ok else "FAIL")
