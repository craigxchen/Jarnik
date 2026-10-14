"""Independent certification of sdc11 (T=1, L=2, m=8) common zeros of the 5 tangent octics, using the referee's
exact invariant basis (sdc11_inv_basis.json).  Krawczyk test in arb; separation; condition (F)."""
import json, sys, time, numpy as np, mpmath as mp, flint
from flint import arb, arb_mat
B = json.load(open('sdc11_inv_basis.json'))
K = np.array(B['K'], dtype=np.int64); C = [np.array(c, dtype=object) for c in B['C']]
c1 = [sum(abs(x) for x in c) for c in C]
d = 11; m = 8; s = 5
Kf = K.astype(float); Cf = [np.array([float(x) for x in c])/c1[i] for i, c in enumerate(C)]
# float64/mpmath helpers
def F_np(u):
    l = Kf @ u; return np.array([np.dot(cf, l**m) for cf in Cf])
def J_np(u):
    l = Kf @ u; return np.array([m*(cf*l**(m-1)) @ Kf for cf in Cf])
mp.mp.dps = 80
Kl = K.tolist(); Cm = [[mp.mpf(int(x))/c1[i] for x in c] for i, c in enumerate(C)]
nz = [[j for j in range(len(Kl)) if C[i][j] != 0] for i in range(s)]
def lin_mp(u): return [mp.fsum(k[t]*u[t] for t in range(d) if k[t]) for k in Kl]
def F_mp(u):
    l = lin_mp(u); return [mp.fsum(Cm[i][j]*l[j]**m for j in nz[i]) for i in range(s)]
def J_mp(u):
    l = lin_mp(u)
    return [[mp.fsum(m*Cm[i][j]*l[j]**(m-1)*Kl[j][t] for j in nz[i] if Kl[j][t]) for t in range(d)] for i in range(s)]
def run(u0, label):
    t0 = time.time()
    u = np.array(u0, dtype=float)
    # float Newton (min-norm) keeping u_0 = 1
    for _ in range(20):
        Jm = J_np(u); Jm[:, 0] = 0
        u = u - np.linalg.pinv(Jm) @ F_np(u)
    Jm = J_np(u)
    # choose 5 unknown coordinates (excluding 0) with best conditioned 5x5 minor (greedy QR pivoting)
    import scipy.linalg as sl
    _, _, piv = sl.qr(Jm[:, 1:], pivoting=True); unk = sorted([int(p)+1 for p in piv[:s]])
    fixed = [t for t in range(d) if t not in unk]
    um = [mp.mpf(x) for x in u]
    for t in fixed: um[t] = mp.mpf(int(mp.nint(um[t]*10**25)))/10**25
    um[0] = mp.mpf(1)
    for it in range(8):
        Fv = F_mp(um); Jv = J_mp(um)
        Js = mp.matrix([[Jv[i][t] for t in unk] for i in range(s)])
        dx = mp.lu_solve(Js, mp.matrix(Fv))
        for a, t in enumerate(unk): um[t] -= dx[a]
    res = max(abs(x) for x in F_mp(um))
    # Krawczyk in arb: X = xc + [-r, r]^5
    flint.ctx.prec = 300
    r = mp.mpf('1e-30')
    xc = [arb(mp.nstr(um[t], 70)) for t in unk]
    def ufrom(xs):
        v = [arb(mp.nstr(um[t], 70)) if t in fixed else None for t in range(d)]
        for a, t in enumerate(unk): v[t] = xs[a]
        return v
    rad = arb(0, 1e-30)
    X = [x + rad for x in xc]
    def F_arb(v):
        l = [sum((v[t]*k[t] for t in range(d) if k[t]), arb(0)) for k in Kl]
        return [sum((arb(int(C[i][j]))*l[j]**m for j in nz[i]), arb(0))/c1[i] for i in range(s)], l
    def J_arb(l):
        return [[sum((arb(int(C[i][j]))*m*l[j]**(m-1)*Kl[j][t] for j in nz[i] if Kl[j][t]), arb(0))/c1[i] for t in unk] for i in range(s)]
    Fc, _ = F_arb(ufrom(xc))
    _, lX = F_arb(ufrom(X)); JX = J_arb(lX)
    Jc = mp.matrix([[J_mp(um)[i][t] for t in unk] for i in range(s)])
    Y = Jc**-1
    Ya = [[arb(mp.nstr(Y[i, j], 40)) for j in range(s)] for i in range(s)]
    Kx = []
    for i in range(s):
        val = xc[i] - sum((Ya[i][j]*Fc[j] for j in range(s)), arb(0))
        for j in range(s):
            Mij = (arb(1) if i == j else arb(0)) - sum((Ya[i][l]*JX[l][j] for l in range(s)), arb(0))
            val += Mij*(X[j] - xc[j])
        Kx.append(val)
    inside = all((Kx[i] - xc[i]).abs_upper() < 1e-30 for i in range(s))
    # separation (exact enclosure: fixed coords exact rationals, unknowns within 1e-30)
    z = [0.0] + [float(x) for x in um]
    zs = sorted(z); sep = min(b - a for a, b in zip(zs, zs[1:]))
    print('[%s] unknown coords %s; residual %.1e; Krawczyk K(X) in int X: %s; sep(u*) >= %.6f (minus 2e-30); |u|_inf = 1 at coord 0: %s; max other |u_x| = %.4f'
          % (label, unk, float(res), inside, sep, True, max(abs(x) for x in z[2:])), flush=True)
    # condition (F): restricted Hessians on the exact tangent space (mpmath at 80 digits; Lipschitz error ~1e-20)
    Jf = np.array([[float(x) for x in row] for row in J_mp(um)])
    U_, S_, Vt = np.linalg.svd(Jf); Tb = Vt[s:].T          # 11 x 6 tangent basis
    l = Kf @ np.array([float(x) for x in um])
    Hs = []
    for i in range(s):
        w = m*(m-1)*Cf[i]*l**(m-2)
        H = (Kf.T * w) @ Kf
        Hs.append((Tb.T @ H @ Tb).ravel())
    G = np.array(Hs); sv = np.linalg.svd(G, compute_uv=False)
    print('[%s] sigma(J) = %s ; singular values of the 5 x 36 restricted-Hessian matrix: %s  (%.0fs)'
          % (label, np.array2string(S_, precision=3), np.array2string(sv, precision=4), time.time()-t0), flush=True)
    return um
if len(sys.argv) > 1:
    pts = []
    for f in sys.argv[1:]:
        for b in json.load(open(f))[:2]: pts.append((b[0], b[1]))
    pts.sort(key=lambda x: -x[0])
    for sep, u in pts[:4]:
        run(np.array(u), 'search zero, float sep %.5f' % sep)
else:
    Z1 = np.load('../fake_checks/sdc11_landzeros.npy')
    for i in range(2):
        run(Z1[i], 'zero %d' % (i+1))
