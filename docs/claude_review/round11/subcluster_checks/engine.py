# Numerical engine (EVIDENCE ONLY) for single and combined base loci of slab-box systems.
# System for (J, base): A = slab-box (|k_x|<=T for x != base, |<k,beta_p|J>|<=w), m = order; V = V(A,m) (full kernel,
# = V_H at these sizes), V^perp = Poly_{<m}|_A.  Residuals:
#   res(delta) = ||(I-P) f_delta||^2/|A|,  f_delta(k) = exp(i<k,delta>)        (base locus Z_A)
#   tau(u)     = ||(I-P) <k,u>^m||^2 / ||<k,u>^m||^2                            (tangent cone)
import numpy as np, itertools, math
from numpy.polynomial import chebyshev as C
from common import slab_box, width, moment_matrix, rank_mod_p

def poly_basis(A, base_idx, m, T, exact_rank=True):
    A = np.asarray(A, dtype=float); n = A.shape[1]
    others = [i for i in range(n) if i != base_idx]
    X = A[:, others]/T
    d = n-1
    mons = [a for t in range(m) for a in itertools.product(range(t+1), repeat=d) if sum(a) == t]
    cols = []
    for al in mons:
        v = np.ones(len(A))
        for j, e in enumerate(al):
            if e: v = v*np.cos(e*np.arccos(np.clip(X[:, j], -1, 1)))
        cols.append(v)
    P = np.stack(cols, axis=1)
    U, s, _ = np.linalg.svd(P, full_matrices=False)
    if exact_rank:
        Mat, _ = moment_matrix([tuple(int(t) for t in k) for k in A], base_idx, m)
        rk = rank_mod_p(Mat)
    else:
        rk = int((s > s[0]*1e-10).sum())
    gap = (s[rk-1]/s[0], s[rk]/s[0] if rk < len(s) else 0.0)
    return U[:, :rk], rk, gap

class System:
    def __init__(self, J, base, cols, T, w, m, exact_rank=True):
        self.J = list(J); self.base = base
        Al = slab_box(self.J, base, cols, T, w)
        self.A = np.array(Al, dtype=float); self.m = m; self.T = T; self.w = w
        bi = self.J.index(base)
        self.Q, self.rk, self.gap = poly_basis(Al, bi, m, T, exact_rank)
        self.dimV = len(Al) - self.rk
        self.h = [width(Al, [2*cols[p][x]-1 for x in self.J]) for p in range(len(cols))]
    def res(self, dJ):
        f = np.exp(1j*(self.A@dJ)); g = f - self.Q@(self.Q.T@f)
        return np.vdot(g, g).real/len(self.A)
    def res_grad(self, dJ):
        f = np.exp(1j*(self.A@dJ)); g = f - self.Q@(self.Q.T@f)
        val = np.vdot(g, g).real/len(self.A)
        grad = (2/len(self.A))*np.real((np.conj(g)*1j*f)@self.A)
        return val, grad
    def tau(self, uJ):
        f = (self.A@uJ)**self.m; nf = np.dot(f, f)
        g = f - self.Q@(self.Q.T@f)
        return np.dot(g, g)/nf

def choose_m(J, base, cols, T, w, omega, min_dimV):
    """smallest m > sum_p omega_p h_A(beta_p) (LP supercritical, leading order); return None if dim V(A,m) < min_dimV"""
    Al = slab_box(list(J), base, cols, T, w)
    h = [width(Al, [2*cols[p][x]-1 for x in J]) for p in range(len(cols))]
    lp = sum(o*hh for o, hh in zip(omega, h))
    m = int(math.floor(lp)) + 1
    Mat, _ = moment_matrix(Al, list(J).index(base), m)
    rk = rank_mod_p(Mat)
    return m, len(Al), len(Al)-rk, lp, max(h)

def int_relation(x, vmax=6, tol=1e-7):
    """smallest (l1) integer v with sum v=0, v != 0, |<v,x>| < tol (x real vector, NOT mod 2pi): exhaustive."""
    n = len(x); best = None
    for v in itertools.product(range(-vmax, vmax+1), repeat=n-1):
        vv = list(v) + [-sum(v)]
        if all(t == 0 for t in vv): continue
        if abs(np.dot(vv, x)) < tol*max(1, np.abs(x).max()):
            if best is None or sum(map(abs, vv)) < sum(map(abs, best)): best = vv
    return best
