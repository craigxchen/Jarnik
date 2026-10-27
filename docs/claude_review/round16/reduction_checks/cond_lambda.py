"""EVIDENCE (float64, SVD): conditioning of the graded form q_A (Corollary 3.2 (Cond)).
For each box A and order m: monomials k^alpha (alpha_x <= 2, |alpha| = m) in the coefficient l2 inner product;
G_top = Gram (in L2(mu_A), probability normalisation) of their top components gr(k^alpha).
Reported: h_A(m) = rank, and the smallest nonzero eigenvalue lambda_min (= min q_A(phi)/|phi|^2 over phi orthogonal to J(A)_m
in coefficient l2). Cube: exact value (2/9)^{floor(m/2)} (2/3)^{m mod 2} for comparison."""
import itertools, sys, time
import numpy as np

def monos(d, deg):
    return [a for a in itertools.product((0, 1, 2), repeat=d) if sum(a) == deg]

def evalm(P, a):
    out = np.ones(len(P))
    for x, ax in enumerate(a):
        if ax:
            out *= P[:, x] ** ax
    return out

def analyse(P, m):
    d = P.shape[1]
    U = len(P)
    low = [a for deg in range(m) for a in monos(d, deg)]
    L = np.stack([evalm(P, a) for a in low], axis=1) / np.sqrt(U)
    Uq, sv, _ = np.linalg.svd(L, full_matrices=False)
    Qb = Uq[:, sv > 1e-9 * sv[0]]
    top = monos(d, m)
    T = np.stack([evalm(P, a) for a in top], axis=1) / np.sqrt(U)
    R = T - Qb @ (Qb.T @ T)
    G = R.T @ R
    ev = np.linalg.eigvalsh(G)
    nz = ev[ev > 1e-10]
    return len(nz), (nz.min() if len(nz) else float('nan')), Qb.shape[1]

def cube(d):
    return np.array(list(itertools.product((-1, 0, 1), repeat=d)), dtype=float)

def onecut(d, S, L):
    C = cube(d)
    return C[np.abs(C[:, S].sum(axis=1)) <= L]

def syl8(L):
    pts = [x for x in itertools.product((0, 1), repeat=3) if x != (0, 0, 0)]  # non-base points; base = 000
    C = cube(7)
    keep = np.ones(len(C), bool)
    for a in itertools.product((0, 1), repeat=3):
        if a == (0, 0, 0):
            continue
        S = [i for i, x in enumerate(pts) if (a[0]*x[0] + a[1]*x[1] + a[2]*x[2]) % 2 == 1]
        keep &= np.abs(C[:, S].sum(axis=1)) <= L
    return C[keep]

out = []
def log(s):
    print(s); out.append(s); sys.stdout.flush()

t0 = time.time()
boxes = [("cube d=6", cube(6)), ("one-cut d=6 S=0..3 L=1", onecut(6, [0, 1, 2, 3], 1)),
         ("one-cut d=6 S=0..3 L=0", onecut(6, [0, 1, 2, 3], 0)),
         ("Syl(8) L=1 (d=7)", syl8(1)), ("M8 = Syl(8) L=2 (d=7)", syl8(2))]
for name, P in boxes:
    for m in range(1, 9):
        h, lam, r0 = analyse(P, m)
        ref = (2/9) ** (m // 2) * (2/3) ** (m % 2)
        log(f"{name:28s} |A|={len(P):5d} m={m}: r0={r0:5d} h_A(m)={h:4d} lambda_min={lam:.3e}  (cube ref {ref:.3e})  [{time.time()-t0:.0f}s]")
open(__file__.replace('.py', '_out.txt'), 'w').write('\n'.join(out) + '\n')
