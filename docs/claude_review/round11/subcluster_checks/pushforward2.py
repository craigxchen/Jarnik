# Sub-cluster systems pushed forward to the coordinates t of a LOW-RANK global profile (subcluster.md Sec. 3-4).
# Profile: points lambda_x in Z^rho; columns y_p = <d_p, lambda> with ranges e_p, weights omega_p (sum 1);
# delta_x = <lambda_x, t>; pi_J(k) = sum_x k_x lambda_x;  slab |<k,beta_p>| = (2/e_p)|<d_p, pi(k)>| <= w.
# Exact (mod p) dimensions of: single pushforwards X_J = pi_*V(A_J,m), combined sum_J X_J, and the universal slab
# space U = V(S_slab, m), S_slab = Lambda cap slab region.  Numerical base loci on T^rho (EVIDENCE).
import sys, itertools, math, numpy as np
from common import *

class Profile:
    def __init__(self, pts, dirs, omega):
        self.pts = [tuple(p) for p in pts]; self.rho = len(pts[0]); self.M = len(pts)
        self.dirs = dirs; self.omega = omega
        ys = [[sum(d[i]*p[i] for i in range(self.rho)) for d in dirs] for p in self.pts]
        self.e = [max(y[j] for y in ys) - min(y[j] for y in ys) for j in range(len(dirs))]
    def norm(self, y):  # y in Z^rho (Lambda coordinates)
        return sum(o*2*abs(sum(d[i]*y[i] for i in range(self.rho)))/e for d, o, e in zip(self.dirs, self.omega, self.e))
    def in_slab(self, y, w):
        return all(2*abs(sum(d[i]*y[i] for i in range(self.rho))) <= w*e for d, e in zip(self.dirs, self.e))

def systems(P, n, T, w, m, bases='all'):
    out = []
    for J in itertools.combinations(range(P.M), n):
        for b in (J if bases == 'all' else [J[-1]]):
            bi = J.index(b); others = [i for i in range(n) if i != bi]
            A = []
            for ks in itertools.product(range(-T, T+1), repeat=n-1):
                k = [0]*n
                for i, v in zip(others, ks): k[i] = v
                k[bi] = -sum(ks)
                y = tuple(sum(k[j]*P.pts[J[j]][i] for j in range(n)) for i in range(P.rho))
                if P.in_slab(y, w): A.append((tuple(k), y))
            hs = []
            for d, e in zip(P.dirs, P.e):
                v = [2*sum(d[i]*a[1][i] for i in range(P.rho))/e for a in A]; hs.append(max(v)-min(v))
            LP = sum(o*h for o, h in zip(P.omega, hs))
            out.append(dict(J=J, base=b, bi=bi, A=A, LP=LP, adm=(m > LP)))
    return out

def poly_rows(Ak, bi, m):
    n = len(Ak[0]); others = [i for i in range(n) if i != bi]
    return np.array([[math.prod(k[o]**a for o, a in zip(others, al)) % P31 for k in Ak] for al in monomials(n-1, m)], dtype=np.int64)

def cheb_poly_float(Ak, bi, m, T):
    A = np.array(Ak, dtype=float); n = A.shape[1]; others = [i for i in range(n) if i != bi]
    X = A[:, others]/max(T, 1)
    cols = []
    for al in monomials(n-1, m):
        v = np.ones(len(A))
        for j, e in enumerate(al):
            if e: v = v*np.cos(e*np.arccos(np.clip(X[:, j], -1, 1)))
        cols.append(v)
    return np.stack(cols, axis=1)

def run(P, n, T, w, m, label, grid=0, bases='all', want_float=True):
    syst = [s for s in systems(P, n, T, w, m, bases) if s['adm']]
    # universal slab support S_slab in Lambda: here Lambda = Z^rho (checked by caller); enumerate a box
    R = max(int(w*max(P.e)) + 2, 2)
    S_slab = [y for y in itertools.product(range(-R, R+1), repeat=P.rho) if P.in_slab(y, w)]
    idx = {y: i for i, y in enumerate(S_slab)}; G = len(S_slab)
    print(f"\n[{label}] rho={P.rho} M={P.M} r={len(P.dirs)} e={P.e} omega={[round(o,3) for o in P.omega]} | n={n} T={T} w={w} m={m}: "
          f"admissible systems {len(syst)}; |S_slab|={G}", flush=True)
    Xs = []; Xf = []; info = []
    if not syst:
        print("  no admissible systems"); return None
    for s in syst:
        Ak = [a[0] for a in s['A']]; Ay = [a[1] for a in s['A']]
        V = nullspace_mod_p(poly_rows(Ak, s['bi'], m))
        Pi = np.zeros((len(Ak), G), dtype=np.int64)
        for j, y in enumerate(Ay): Pi[j, idx[y]] = 1
        X = (V @ Pi) % P31 if len(V) else np.zeros((0, G), dtype=np.int64)
        rX = rank_rows_mod_p(X) if len(X) else 0
        info.append((len(Ak), len(V), rX, len(set(Ay)), s['LP']))
        if rX:
            Xs.append(X)
            if want_float:
                Pf = cheb_poly_float(Ak, s['bi'], m, T)
                U_, sv, _ = np.linalg.svd(Pf, full_matrices=True)
                rk = len(Ak) - len(V)
                Vf = U_[:, rk:]                       # orthonormal basis of V(A,m) (float)
                Xfl = Pi.astype(float).T @ Vf          # G x dimV
                Uq, sq, _ = np.linalg.svd(Xfl, full_matrices=False)
                Xf.append((Uq[:, :rX], sq[rX-1]/sq[0], sq[rX]/sq[0] if rX < len(sq) else 0.0))
    dim_comb = rank_rows_mod_p(np.concatenate(Xs)) if Xs else 0
    Pbox = np.array([[math.prod(g[i]**al[i] for i in range(P.rho)) % P31 for g in S_slab] for al in monomials(P.rho, m)], dtype=np.int64)
    dim_U = G - rank_rows_mod_p(Pbox)
    singles = [i[2] for i in info]
    print(f"  per system: |A| {min(i[0] for i in info)}..{max(i[0] for i in info)}, dim V(A,m) {min(i[1] for i in info)}..{max(i[1] for i in info)}, "
          f"dim pi_*V {min(singles)}..{max(singles)} (zero for {sum(1 for x in singles if x==0)}), LP {min(i[4] for i in info):.2f}..{max(i[4] for i in info):.2f}")
    print(f"  dim combined = {dim_comb};  dim universal U = V(S_slab,m) = {dim_U}", flush=True)
    return dict(syst=syst, Xf=Xf, S=S_slab, dim_comb=dim_comb, dim_U=dim_U, singles=singles, info=info)

def float_space(Xf_list):
    if not Xf_list: return None
    B = np.concatenate([x[0] for x in Xf_list], axis=1)
    U_, s, _ = np.linalg.svd(B, full_matrices=False)
    return U_, s

def residual_fn(Q, S):
    Sa = np.array(S, dtype=float)
    def f(t):
        chi = np.exp(1j*(Sa @ t)); c = Q.T @ chi
        return np.vdot(c, c).real/len(Sa)
    return f

if __name__ == "__main__":
    # rank-2 checkerboard profile (5 points, D2 lattice): lambda in {(0,0),(1,0),(1,1),(1,2),(2,2)}, columns (-2,1),(0,-1)
    P = Profile([(0,0),(1,0),(1,1),(1,2),(2,2)], [(-2,1),(0,-1)], [0.5,0.5])
    print("pair norms:", sorted(round(P.norm(np.subtract(a,b)),3) for a,b in itertools.combinations(P.pts,2)))
    for (n, T, w, m) in [(4,2,2,5),(4,3,2,5),(4,3,2,6),(5,2,2,5),(5,2,2,6),(4,3,3,7),(4,3,3,8),(5,2,3,7),(5,3,3,7)]:
        run(P, n, T, w, m, "checkerboard rho=2", want_float=False)
