"""Referee checks on the record clusters of:
 (i)  Thm A of p3lower (W >= 2 lambda_max(D) - 2(M-1) c_C^+) versus Prop 2.2 of L1L2 (lambda_max(D) <= max_x D_x
      <= W/2 + sqrt(W Q) + (M-1) c_C), with exact D at moduli q <= X (q odd, q not dividing N');
 (ii) Prop 5.1(b) area identity 4(2 Area)^2 = Norm(g_xyw)^2 |b_xy|^2 |b_yw|^2 |b_xw|^2 (exact);
 (iii) Prop 5.2 (short relation sum c_x z_x = 0, sum c_x = 0; Minkowski bound) via sympy LLL, and Lemma 5.3
      (threshold cuts meet supp c in 0 or >= 2 points on each side, for columns p_j not dividing any c_x).
"""
import itertools
from math import log, sqrt
import numpy as np
from refgauss import *
from records import profile, primes_upto
from sympy import Matrix
from sympy.polys.matrices import DomainMatrix
from sympy import ZZ

targets = [("9pt", 1176852625, 9, 8.0), ("10pt", 1176852625, 10, 8.3), ("7pt", 235370525, 7, 5.3),
           ("8pt golden", 537606725, 8, 5.9)]

for name, N, M, Cm in targets:
    ws = sorted(windows(N, M, Cm), key=lambda t: t[1])
    cl, C = ws[0]
    u, G = normalise(cl)
    Np, fac, pis, rows, ks = profile(u)
    W = log(Np)
    cC = 2 * log(C) - log(2)
    # (i) exact D with X = 10^4
    X = 10 ** 4
    qs = [q for q in primes_upto(X) if q > 2 and Np % q]
    D = np.zeros((M, M)); mu = np.zeros((M, M))
    for x, y in itertools.combinations(range(M), 2):
        g = ggcd(u[x], u[y]); beta = gexactdiv(gsub(u[x], u[y]), g)
        # m_xy: largest odd q-power product with u_x = u_y mod q^a, q <= X, q not | N'
        m = 1
        for q in qs:
            Q = q
            while (u[x][0] - u[y][0]) % Q == 0 and (u[x][1] - u[y][1]) % Q == 0:
                m *= q; Q *= q
        D[x, y] = D[y, x] = 2 * log(m)
        mu[x, y] = mu[y, x] = log(Np) - log(gnorm(g)) - W / 2
    lam = np.linalg.eigvalsh(D)[-1]
    Dx = D.sum(1).max()
    # cut slack Q = sum_cuts w (S - M/2)^2
    Qs = 0.0
    for p, e in fac.items():
        for t in range(1, e + 1):
            S = sum(1 for r in rows if r[p] < t)
            Qs += log(p) * (S - M / 2) ** 2
    print(f"{name}: W={W:.3f} c_C={cC:.3f} lambda_max(D)={lam:.3f} max_x D_x={Dx:.3f} | ThmA ceiling W/2+(M-1)c_C^+="
          f"{W/2+(M-1)*max(cC,0):.3f} | Prop2.2 first line W/2+sqrt(WQ)+(M-1)c_C={W/2+sqrt(W*Qs)+(M-1)*cC:.3f}  (Q={Qs:.3f})")
    assert W >= 2 * lam - 2 * (M - 1) * max(cC, 0) - 1e-9
    # (ii) area identity, all triples
    nt = 0
    for x, y, w in itertools.combinations(range(M), 3):
        a, b, c = u[x], u[y], u[w]
        area2 = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))  # = 2 Area
        g3 = ggcd(ggcd(a, b), c)
        bxy = gnorm(gsub(a, b)) // gnorm(ggcd(a, b))
        byw = gnorm(gsub(b, c)) // gnorm(ggcd(b, c))
        bxw = gnorm(gsub(a, c)) // gnorm(ggcd(a, c))
        assert 4 * area2 ** 2 == gnorm(g3) ** 2 * bxy * byw * bxw
        nt += 1
    # (iii) relation lattice: kernel of c -> (sum c, sum c Re z, sum c Im z)
    A = Matrix([[1] * M, [z[0] for z in u], [z[1] for z in u]])
    ker = A.nullspace()
    # integer basis: clear denominators then saturate via HNF-free approach: use LLL on a big-weight embedding
    BIG = 10 ** 30
    rowsL = []
    for i in range(M):
        e_ = [0] * M; e_[i] = 1
        rowsL.append(e_ + [BIG * 1, BIG * u[i][0], BIG * u[i][1]])
    dm = DomainMatrix([[ZZ(v) for v in r] for r in rowsL], (M, M + 3), ZZ)
    red = dm.lll().to_Matrix()
    rels = []
    for i in range(M):
        r = list(red.row(i))
        if all(v == 0 for v in r[M:]):
            rels.append(r[:M])
    assert len(rels) == M - 3, len(rels)
    covol = sqrt(abs(Matrix(rels) * Matrix(rels).T).det())
    shortest = min(sqrt(sum(v * v for v in r)) for r in rels)
    mink = sqrt(M - 3) * (M ** 1.5 * C ** 2 * sqrt(Np)) ** (1 / (M - 3))
    print(f"   area identity on {nt} triples: OK.  relation lattice rank {len(rels)}, covol={covol:.4g} "
          f"(bound M^1.5 C^2 N^.5 = {M**1.5*C**2*sqrt(Np):.4g}); shortest LLL |c|_2={shortest:.3f}, Minkowski bound={mink:.3f}")
    # Lemma 5.3 on all LLL relations
    bad = 0; chk = 0
    for c in rels:
        supp = [x for x in range(M) if c[x] != 0]
        for p, e in fac.items():
            if any(c[x] % p == 0 for x in supp):
                continue
            for t in range(1, e + 1):
                lo = sum(1 for x in supp if rows[x][p] < t); hi = len(supp) - lo
                chk += 1
                if lo == 1 or hi == 1:
                    bad += 1
    print(f"   Lemma 5.3: {chk} (relation, column, threshold) checks, violations {bad}")
