# kill_search: Krawczyk certificate for a rank-4, 5-point cluster on the vertical line q_k = i y_k, D=10, c=0.
import numpy as np, itertools, mpmath as mp, networkx as nx
mp.mp.dps = 60; iv = mp.iv; iv.dps = 60
sols = np.load('vert_sols.npy'); D = 10; L = 4
subs = np.array(list(itertools.product([0,1], repeat=D)), float)
for idx in [1, 2]:
    y = sols[idx]
    R = np.sqrt((subs @ y)**2 + (subs @ y**3)**2)
    zs = [s for s, rv in zip(subs, R) if rv < 1e-9]
    basis = [s for s in zs if s.sum() > 0][:4]
    assert np.linalg.matrix_rank(np.array(basis)) == 4
    j0 = int(np.argmax(np.abs(y)))  # coordinate to fix
    fixval = mp.mpf(round(y[j0]*10**8))/10**8
    def f(u):
        out = []
        for b in basis:
            out.append(mp.fsum(u[k] for k in range(D) if b[k]))
            out.append(mp.fsum(u[k]**3 for k in range(D) if b[k]))
        out.append(mp.fsum(x**2 for x in u) - D)
        out.append(u[j0] - fixval)
        return out
    def J(u):
        Jm = mp.matrix(D, D)
        for i, b in enumerate(basis):
            for k in range(D):
                if b[k]: Jm[2*i, k] = 1; Jm[2*i+1, k] = 3*u[k]**2
        for k in range(D): Jm[8, k] = 2*u[k]
        Jm[9, j0] = 1
        return Jm
    u = [mp.mpf(float(v)) for v in y]
    for it in range(60):
        fv = mp.matrix(f(u))
        if mp.norm(fv) < mp.mpf(10)**-55: break
        dx = mp.lu_solve(J(u), fv); u = [u[k]-dx[k] for k in range(D)]
    rad = mp.mpf(10)**-25
    Y = mp.inverse(J(u))
    X = [iv.mpf([x-rad, x+rad]) for x in u]
    JX = [[iv.mpf(0)]*D for _ in range(D)]
    for i, b in enumerate(basis):
        for k in range(D):
            if b[k]: JX[2*i][k] = iv.mpf(1); JX[2*i+1][k] = 3*X[k]**2
    for k in range(D): JX[8][k] = 2*X[k]
    JX[9][j0] = iv.mpf(1)
    uc = [iv.mpf(x) for x in u]
    fc = []
    for b in basis:
        fc.append(sum((uc[k] for k in range(D) if b[k]), iv.mpf(0)))
        fc.append(sum((uc[k]**3 for k in range(D) if b[k]), iv.mpf(0)))
    fc.append(sum((x**2 for x in uc), iv.mpf(0)) - D); fc.append(uc[j0] - iv.mpf(fixval))
    ok = True
    for i in range(D):
        val = uc[i] - sum((iv.mpf(Y[i, l])*fc[l] for l in range(D)), iv.mpf(0))
        for col in range(D):
            coef = (1 if i == col else 0) - sum((iv.mpf(Y[i, l])*JX[l][col] for l in range(D)), iv.mpf(0))
            val += coef*(X[col] - uc[col])
        if not (val.a > X[i].a and val.b < X[i].b): ok = False
    yy = np.array([float(x) for x in u])
    ay = np.sort(np.abs(yy))
    # full model check: q = i*y, all moments m=1..L, cluster through Q and max cluster over all points of N
    q = 1j*yy
    A = np.array([np.imag(q**m) for m in range(1, L+1)])
    res = np.linalg.norm(subs @ A.T, axis=1)
    nQ = int((res < 1e-9).sum())
    W = np.array(list(itertools.product([-1,0,1], repeat=D)), float)
    rr = np.linalg.norm(W @ A.T, axis=1); ker = W[rr < 1e-9].astype(int)
    G = nx.Graph(); G.add_nodes_from(range(2**D))
    for i, s in enumerate(subs.astype(int)):
        for kv in ker:
            g = s - kv
            if np.all((g == 0) | (g == 1)):
                jj = int(''.join(map(str, g)), 2)
                if jj > i: G.add_edge(i, jj)
    cl, _ = nx.max_weight_clique(G, weight=None)
    print("sol", idx, "KRAWCZYK", "PASS" if ok else "FAIL", "min|y|", ay[0], "min gap of |y|", np.min(np.diff(ay)),
          "cluster through Q", nQ, "rank", np.linalg.matrix_rank(subs[res < 1e-9]), "min nonmember residual", np.min(res[res >= 1e-9]),
          "max cluster on N", len(cl))
    print("  y =", [mp.nstr(x, 20) for x in u])
    print("  basis sets:", [list(np.where(b)[0]) for b in basis])
