"""NUMERICAL: Klein ansatz in t-coordinates (see klein.py for the ansatz).
Unknowns: t[x][i] (x in {0,1,inf}, i=1..4): representative of orbit O_i^x where a=b=x;
t_cd, t_abc, t_abd; kappa_A, kappa_B, kappa_C, kappa_D; omega (barrier).
"""
import numpy as np, sys

def esym(r):
    # elementary symmetric polys of roots r (monic poly coefficients without leading 1)
    return np.poly(r)[1:]

def unpack(u):
    t = u[0:12].reshape(3, 4)
    tcd, tabc, tabd = u[12], u[13], u[14]
    kA, kB, kC, kD = u[15], u[16], u[17], u[18]
    om = u[19]
    return t, tcd, tabc, tabd, kA, kB, kC, kD, om

def funcs(u):
    t, tcd, tabc, tabd, kA, kB, kC, kD, om = unpack(u)
    v = t + 1/t; y = t**2; z = y + 1/y
    PA = [np.poly(v[x]) for x in range(3)]
    PB = [np.poly(y[x]) for x in range(3)]
    QC = [np.poly(z[x, [0, 1]]) for x in range(3)]
    QD = [np.poly(z[x, [1, 2]]) for x in range(3)]
    A = lambda w: kA*np.polyval(PA[0], w)/np.polyval(PA[2], w)
    B = lambda w: kB*np.polyval(PB[0], w)/np.polyval(PB[2], w)
    C = lambda w: kC*np.polyval(QC[0], w)/np.polyval(QC[2], w)
    D = lambda w: kD*np.polyval(QD[0], w)/np.polyval(QD[2], w)
    return t, v, y, z, PA, PB, QC, QD, A, B, C, D

def residual(u):
    t, tcd, tabc, tabd, kA, kB, kC, kD, om = unpack(u)
    t_, v, y, z, PA, PB, QC, QD, A, B, C, D = funcs(u)
    R = []
    for P, k in [(PA, kA), (PB, kB), (QC, kC), (QD, kD)]:
        dep = k*P[0] - P[2] - (k-1)*P[1]
        R += list(dep[1:])
    def vyz(s):
        return s + 1/s, s**2, s**2 + s**-2
    vv, yy, zz = vyz(tcd)
    g = C(zz)
    R += [D(zz) - g, A(vv) - g, B(yy) - g]
    vv, yy, zz = vyz(tabc)
    g = C(zz)
    R += [A(vv) - g, B(yy) - g]
    vv, yy, zz = vyz(tabd)
    g = D(zz)
    R += [A(vv) - g, B(yy) - g]
    # barrier: all 15 z-values distinct and away from {2,-2}
    zs = list(z.ravel()) + [tcd**2 + tcd**-2, tabc**2 + tabc**-2, tabd**2 + tabd**-2]
    s = 0
    for i in range(len(zs)):
        s = s + np.log(zs[i] - 2) + np.log(zs[i] + 2)
        for j in range(i+1, len(zs)):
            s = s + np.log(zs[i] - zs[j])
    r = s + om
    R.append(r.real + 1j*((r.imag + np.pi) % (2*np.pi) - np.pi))
    return np.array(R)

def newton(u, iters=80, tol=1e-12):
    n = len(u)
    for it in range(iters):
        F = residual(u)
        nF = np.linalg.norm(F)
        if not np.isfinite(nF): return None, it
        if nF < tol: return u, it
        J = np.zeros((n, n), complex); h = 1e-7
        for k in range(n):
            du = np.zeros(n, complex); du[k] = h
            J[:, k] = (residual(u+du) - residual(u-du))/(2*h)
        try: du = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError: return None, it
        st = 1.0
        while st > 1e-9:
            u1 = u + st*du; F1 = residual(u1)
            if np.all(np.isfinite(F1)) and np.linalg.norm(F1) < nF*(1-0.05*st): break
            st /= 2
        if st <= 1e-9: return None, it
        u = u1
    return None, iters

if __name__ == '__main__':
    nst = int(sys.argv[1]); start = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    for seed in range(start, start+nst):
        rng = np.random.default_rng(seed)
        u = rng.normal(size=20) + 1j*rng.normal(size=20)
        u, it = newton(u)
        if u is None:
            print(seed, 'fail', it, flush=True); continue
        np.save('k3sol_%d.npy' % seed, u)
        print(seed, 'converged', it, flush=True)
