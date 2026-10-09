"""NUMERICAL: balanced rational curves in Mbar_{0,7} with one single-label degree-two forgetful map,
built as (balanced rational 6-curve in y) + (double cover y = (y+ tau^2 - y-)/(tau^2 - 1)) + (7th point
d(tau) of degree 8 meeting each of the 25 six-point contact clusters at exactly one of the two lifts).
The 'one of two lifts' condition is imposed symmetrically: (d(tau)-v)(d(-tau)-v) = 0, a polynomial in
s = tau^2 (no choice of pattern).  Expected dimension 2 (Section 5 of L7.md)."""
import numpy as np, itertools, sys
from gen6 import residual6, cluster_values, mover_polys, KEYS, KIDX, MOV, ANC, SUBS, QSUB

FIX = [0, 1, 2]          # indices of y's frozen (PGL_2 normalization)

class Sys:
    def __init__(self, v0, rng):
        self.v0 = v0.copy()
        self.free = [k for k in range(28) if k not in FIX]
        self.L = rng.normal(size=18) + 1j*rng.normal(size=18)
    def split(self, u):
        v = self.v0.copy(); v[self.free] = u[:25]
        yp, ym = u[25], u[26]
        a = u[27:45]
        return v, yp, ym, a
    def residual(self, u):
        v, yp, ym, a = self.split(u)
        R = list(residual6(v))
        N = a[:9]; D = a[9:]          # coefficients highest degree first, degree 8
        # even/odd parts in s = tau^2 : N(tau) = Ne(s) + tau No(s)
        def eo(P):
            c = P[::-1]               # c[k] = coeff of tau^k
            Ne = c[0::2][::-1]; No = c[1::2][::-1]
            return Ne, No
        Ne, No = eo(N); De, Do = eo(D)
        y = v[:25]; vals = cluster_values(v)
        for j in range(25):
            s = (y[j] - ym)/(y[j] - yp)
            if vals[j] is None:
                e = np.polyval(De, s)**2 - s*np.polyval(Do, s)**2
            else:
                vv = vals[j]
                e = (np.polyval(Ne, s) - vv*np.polyval(De, s))**2 - s*(np.polyval(No, s) - vv*np.polyval(Do, s))**2
            R.append(e/(1 + abs(s))**8)
        R.append(np.dot(self.L, a) - 1)
        return np.array(R)

def gauss_newton(S, u, iters=150, tol=1e-12):
    n = len(u)
    for it in range(iters):
        F = S.residual(u); nF = np.linalg.norm(F)
        if not np.isfinite(nF): return None, it
        if nF < tol: return u, it
        J = np.zeros((len(F), n), complex); h = 1e-7
        for k in range(n):
            du = np.zeros(n, complex); du[k] = h
            J[:, k] = (S.residual(u+du) - S.residual(u-du))/(2*h)
        du = np.linalg.lstsq(J, -F, rcond=None)[0]
        st = 1.0
        while st > 1e-10:
            F1 = S.residual(u + st*du)
            if np.all(np.isfinite(F1)) and np.linalg.norm(F1) < nF*(1-0.05*st): break
            st /= 2
        if st <= 1e-10: return None, it
        u = u + st*du
    return None, iters

if __name__ == '__main__':
    nst = int(sys.argv[1]); start = int(sys.argv[2])
    v0 = np.load('gen6_seed.npy')
    for seed in range(start, start+nst):
        rng = np.random.default_rng(seed)
        S = Sys(v0, rng)
        u = np.concatenate([v0[S.free] + 0.05*(rng.normal(size=25) + 1j*rng.normal(size=25)),
                            rng.normal(size=2) + 1j*rng.normal(size=2),
                            rng.normal(size=18) + 1j*rng.normal(size=18)])
        r, it = gauss_newton(S, u)
        print(seed, 'OK' if r is not None else 'fail', it, flush=True)
        if r is not None:
            np.save('lift7_sol_%d.npy' % seed, r); np.save('lift7_L_%d.npy' % seed, S.L)
