# Independent certifier (kill_search). Written from scratch; does not import lane code.
# Usage: kcert.py file.npy row c
import sys, itertools, numpy as np, mpmath as mp
from fractions import Fraction
import sympy as sp
import networkx as nx

mp.mp.dps = 60
import os
DEBUG=os.environ.get('DEBUG')

def Lof(D, c):
    K = (D + 2*c)//2  # floor(D/2 + c) for integer c
    return D - 1 - K

def moments_np(q, L):
    return np.array([np.imag(q**m) for m in range(1, L+1)])  # L x D

def load(fn, row):
    a = np.load(fn)
    if a.ndim == 1: return a
    return a[row]

def main():
    fn, row, c = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    q = load(fn, row).astype(complex); D = len(q); L = Lof(D, c)
    print("file", fn, "row", row, "D", D, "c", c, "L", L)
    A = moments_np(q, L)
    subs = list(itertools.product([0,1], repeat=D))
    V = np.array(subs, dtype=float)
    res = np.linalg.norm(V @ A.T, axis=1)
    zs = [i for i in range(len(subs)) if res[i] < 1e-8]
    print("numerical zero-sum sets:", len(zs), " next smallest residual:", np.sort(res[res >= 1e-8])[:3])
    Z = V[zs]
    r = np.linalg.matrix_rank(Z, tol=1e-9)
    print("rank", r)
    # basis of zero-sum sets
    basis = []
    for i in zs:
        if np.linalg.matrix_rank(np.array(basis + [V[i]]), tol=1e-9) > len(basis):
            basis.append(V[i])
    B = np.array(basis)
    # exact check: every numerical zero-sum set lies in the rational span of basis
    Bs = sp.Matrix(B.astype(int))
    for i in zs:
        M2 = Bs.col_join(sp.Matrix([list(V[i].astype(int))]))
        assert M2.rank() == Bs.rank()
    # unknowns u = (x_1..x_D, y_1..y_D); equations sum_{F in basis} Im q^m
    nunk = 2*D; neq = r*L
    d = nunk - neq
    print("free dims (incl. 2 affine):", d)
    # jacobian numerical
    def F_mp(u):
        qs = [mp.mpc(u[k], u[D+k]) for k in range(D)]
        out = []
        for b in basis:
            for m in range(1, L+1):
                out.append(mp.fsum(mp.im(qs[k]**m) for k in range(D) if b[k] > 0.5))
        return out
    def J_np(qv):
        J = np.zeros((neq, nunk))
        row_ = 0
        for b in basis:
            for m in range(1, L+1):
                der = m*qv**(m-1)
                for k in range(D):
                    if b[k] > 0.5:
                        J[row_, k] = np.imag(der[k]); J[row_, D+k] = np.real(der[k])
                row_ += 1
        return J
    J = J_np(q)
    # choose d coordinates to fix: complement of best-conditioned column subset (QR pivoting)
    import scipy.linalg as sl
    _, _, piv = sl.qr(J, pivoting=True)
    keep = sorted(piv[:neq]); fix = sorted(piv[neq:])
    print("fixed coords", fix, "cond of square part", np.linalg.cond(J[:, keep]))
    u0 = [mp.mpf(float(np.real(q[k]))) for k in range(D)] + [mp.mpf(float(np.imag(q[k]))) for k in range(D)]
    for j in fix:  # fix to rationals with 8 digits
        u0[j] = mp.mpf(round(float(u0[j])*10**8)) / 10**8
    # Newton on keep coords
    u = list(u0)
    def Jmp(u):
        qs = [mp.mpc(u[k], u[D+k]) for k in range(D)]
        Jm = mp.matrix(neq, neq); row_ = 0
        for b in basis:
            for m in range(1, L+1):
                for col, j in enumerate(keep):
                    k = j % D
                    if b[k] > 0.5:
                        der = m*qs[k]**(m-1)
                        Jm[row_, col] = mp.im(der) if j < D else mp.re(der)
                row_ += 1
        return Jm
    for it in range(40):
        f = mp.matrix(F_mp(u))
        if DEBUG: print(it, mp.nstr(mp.norm(f),3))
        if mp.norm(f) < mp.mpf(10)**(-55): break
        dx = mp.lu_solve(Jmp(u), f)
        for col, j in enumerate(keep): u[j] -= dx[col]
    print("Newton residual", mp.nstr(mp.norm(mp.matrix(F_mp(u))), 5), "move", mp.nstr(max(abs(u[j]-u0[j]) for j in keep), 5))
    # Krawczyk with mpmath.iv
    iv = mp.iv; iv.dps = 60
    rad = mp.mpf(10)**(-25)
    Y = mp.inverse(Jmp(u))
    X = [iv.mpf([u[j]-rad, u[j]+rad]) if j in keep else iv.mpf(u[j]) for j in range(nunk)]
    def ipow(x, y, m):  # (x+iy)^m as interval pair
        re, im = iv.mpf(1), iv.mpf(0)
        for _ in range(m):
            re, im = re*x - im*y, re*y + im*x
        return re, im
    # interval Jacobian on box
    JX = [[iv.mpf(0)]*neq for _ in range(neq)]
    row_ = 0
    for b in basis:
        for m in range(1, L+1):
            for col, j in enumerate(keep):
                k = j % D
                if b[k] > 0.5:
                    re, im = ipow(X[k], X[D+k], m-1)
                    JX[row_][col] = m*im if j < D else m*re
            row_ += 1
    # f at centre as tight interval
    Xc = [iv.mpf(u[j]) for j in range(nunk)]
    fc = []
    for b in basis:
        for m in range(1, L+1):
            s = iv.mpf(0)
            for k in range(D):
                if b[k] > 0.5:
                    s += ipow(Xc[k], Xc[D+k], m)[1]
            fc.append(s)
    ok = True; maxw = 0
    for i in range(neq):
        val = iv.mpf(u[keep[i]]) - sum(iv.mpf(Y[i, l])*fc[l] for l in range(neq))
        for col in range(neq):
            coef = (1 if i == col else 0) - sum(iv.mpf(Y[i, l])*JX[l][col] for l in range(neq))
            val += coef*(X[keep[col]] - iv.mpf(u[keep[col]]))
        lo, hi = X[keep[i]].a, X[keep[i]].b
        if not (val.a > lo and val.b < hi): ok = False
        maxw = max(maxw, float(mp.mpf((val.b - val.a).a)/(2*rad)) if hasattr(val.b-val.a,'a') else 0)
    print("KRAWCZYK", "PASS" if ok else "FAIL", "contraction ratio", maxw)
    qs = np.array([complex(float(u[k]), float(u[D+k])) for k in range(D)])
    allr = np.concatenate([qs, qs.conj()])
    sep = min(abs(allr[i]-allr[j]) for i in range(2*D) for j in range(i))
    print("min|Im q|", np.min(abs(qs.imag)), "separation of 2D roots", sep)
    # count cluster through Q at certified point
    A = moments_np(qs, L); res = np.linalg.norm(V @ A.T, axis=1)
    zs2 = [i for i in range(len(subs)) if res[i] < 1e-8]
    print("cluster through Q at certified point:", len(zs2), "min nonmember residual", np.min(res[res >= 1e-8]))
    # full max clique on all 2^D points: kernel vectors in {0,+-1}^D
    W = np.array(list(itertools.product([-1,0,1], repeat=D)), dtype=np.int8)
    rr = np.linalg.norm(W.astype(float) @ A.T, axis=1)
    ker = W[rr < 1e-8]
    print("kernel vectors in {0,+-1}^D (incl 0):", len(ker), "min nonzero residual", np.min(rr[rr >= 1e-8]))
    # graph: F~G iff 1_F-1_G in ker
    kerset = set(map(tuple, ker.tolist()))
    G = nx.Graph(); G.add_nodes_from(range(len(subs)))
    for i, s in enumerate(subs):
        for a in ker:
            g = np.array(s) - a
            if np.all((g == 0) | (g == 1)):
                j = int(''.join(map(str, g.tolist())), 2)
                if j > i: G.add_edge(i, j)
    w, _ = nx.max_weight_clique(G, weight=None)
    print("max cluster over ALL points of N:", len(w))

main()
