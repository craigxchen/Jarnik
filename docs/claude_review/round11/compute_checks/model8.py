"""Generic faithful-model builder at n = 8 (T = 1): profile by name, exact rho, LP bound on zeta, exact basis.
Profiles: had8 (Sylvester), sdc7 (skew double core at q = 7, the fake profile of round5/exact.md Thm 4.1),
pal7 (round-5 Paley-I fake profile at q = 7, distinct patterns), rand<seed> (random squarefree sign profile).
Usage: model8.py <profile> T L [m]   (m defaults to rho(A); the model is faithful iff LPbound < m)."""
import sys, json, random, time
import numpy as np
import flint
from core import build_A, moment_matrix, width, siegel_H
from models import canon, fullrank
from lpzeta import lp_bound
from bigger import hadamard8

def legendre(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)

def paley_I(q):
    idx = ['inf'] + list(range(q)); M = q + 1
    H = np.zeros((M, M), dtype=int)
    for a, x in enumerate(idx):
        for b, y in enumerate(idx):
            if x == 'inf': H[a, b] = 1
            elif y == 'inf': H[a, b] = -1
            else: H[a, b] = (1 if x == y else 0) + legendre(y - x, q)
    return H

def sdc(q):
    K = paley_I(q); M = q + 1
    if not ((K + K.T) == 2 * np.eye(M, dtype=int)).all():
        K = K.T
    Dg = np.diag(K[:, 0]); K = Dg @ K @ Dg
    pats = set()
    for a in range(1, M):
        S = [x for x in range(M) if K[x, a] == -1]
        Sp = sorted(set(S) ^ {0, a})
        pats.add(canon(M, tuple(S))); pats.add(canon(M, tuple(Sp)))
    return tuple(sorted(pats))

def pal(q):
    H = paley_I(q); M = q + 1
    Hn = H.copy(); Hn[1:, :] *= -1
    lab = {0: 1}
    for i in range(q): lab[1 + i] = 1 + i
    pats = set()
    for a in range(1, M):
        S = tuple(x for x in range(M) if Hn[x, a] == 1)
        if 0 < len(S) < M: pats.add(canon(M, S))
    for x in range(M):
        c = Hn[:, lab[x]].copy(); c[x] *= -1
        S = tuple(y for y in range(M) if c[y] == 1)
        if 0 < len(S) < M: pats.add(canon(M, S))
    return tuple(sorted(pats))

def profile(name):
    if name == 'had8': return hadamard8()
    if name == 'sdc7': return sdc(7)
    if name == 'pal7': return pal(7)
    if name.startswith('rand'):
        rng = random.Random(int(name[4:]))
        while True:
            r = rng.randint(8, 24); pats = set()
            for _ in range(r):
                S = tuple(x for x in range(8) if rng.random() < 0.5)
                if 0 < len(S) < 8: pats.add(canon(8, S))
            p = tuple(sorted(pats))
            if fullrank(8, p): return p

if __name__ == "__main__":
    name = sys.argv[1]; T = int(sys.argv[2]); L = int(sys.argv[3])
    prof = profile(name); n = 8
    A = build_A(n, T, [(S, L) for S in prof])
    rngs = [width(A, S) for S in prof]
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    pr("%s T=%d L=%d |A|=%d #patterns=%d ranges=%s" % (name, T, L, len(A), len(prof), rngs))
    U = len(A); j = 0
    t0 = time.time()
    while moment_matrix(A, j + 1).rank() < U:
        j += 1
    rho = j
    lb, _ = lp_bound(n, T, prof, L)
    pr("rho(A) = %d (exact), LP bound zeta <= %s (%.1fs); faithful m: %s" % (rho, lb, time.time() - t0,
       [mm for mm in range(int(lb) + 1, rho + 1)]))
    m = int(sys.argv[4]) if len(sys.argv) > 4 else rho
    if not (lb < m <= rho) and "force" not in sys.argv:
        pr("m = %d not faithful; stop" % m); sys.exit(0)
    M = moment_matrix(A, m)
    X, nul = M.nullspace()
    B = flint.fmpz_mat([[X[i, jj] for i in range(U)] for jj in range(nul)]).lll()
    basis = [[int(B[i, jj]) for jj in range(U)] for i in range(B.nrows())]
    for c in basis:
        v = M * flint.fmpz_mat([[x] for x in c]); assert all(v[i, 0] == 0 for i in range(v.nrows()))
    sup = max(abs(x) for c in basis for x in c)
    pr("m=%d: dim V(A,m) = %d, LLL basis sup-norm %d, H = %.3g, V_H = V: %s; all verified exactly" % (
        m, len(basis), sup, siegel_H(A, T, m), sup <= siegel_H(A, T, m)))
    pr("max h_A(beta_p) = %d vs m = %d: (ii) %s, (ii') %s" % (2 * max(rngs), m, 2 * max(rngs) < m, max(rngs) < m))
    tag = "%s_T%d_L%d_m%d" % (name, T, L, m)
    json.dump(dict(A=A, basis=basis, T=T, L=L, m=m, prof=prof), open(tag + "_basis.json", "w"))
    open(tag + "_out.txt", "w").write("\n".join(out) + "\n")
