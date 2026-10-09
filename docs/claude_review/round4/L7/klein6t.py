"""NUMERICAL calibration: k=6 Klein ansatz in t-coordinates with the same barrier as klein3.py.
Unknowns: t[x][0] (orbit with c=x, a=b=x), t[x][1] (orbit with a=b=x, c!=x), t_m (a=b=c), kA, kB, omega;
C is the Moebius map determined by C(z(t[x][0])) = x.  Three random fixed values: t[x][0] (x=0,1,inf)."""
import numpy as np, sys

def residual(u, fixed):
    t = np.zeros((3, 2), complex); t[:, 0] = fixed; t[:, 1] = u[0:3]
    tm, kA, kB, om = u[3], u[4], u[5], u[6]
    v = t + 1/t; y = t**2; z = y + 1/y
    PA = [np.poly(v[x]) for x in range(3)]; PB = [np.poly(y[x]) for x in range(3)]
    R = []
    for P, k in [(PA, kA), (PB, kB)]:
        dep = k*P[0] - P[2] - (k-1)*P[1]; R += list(dep[1:])
    zc = z[:, 0]
    C = lambda w: (w - zc[0])*(zc[1] - zc[2])/((w - zc[2])*(zc[1] - zc[0]))
    A = lambda w: kA*np.polyval(PA[0], w)/np.polyval(PA[2], w)
    B = lambda w: kB*np.polyval(PB[0], w)/np.polyval(PB[2], w)
    vv, yy, zz = tm + 1/tm, tm**2, tm**2 + tm**-2
    R += [A(vv) - C(zz), B(yy) - C(zz)]
    zs = list(z.ravel()) + [zz]
    sm = 0
    for i in range(len(zs)):
        sm = sm + np.log(zs[i] - 2) + np.log(zs[i] + 2)
        for j in range(i+1, len(zs)): sm = sm + np.log(zs[i] - zs[j])
    r = sm + om
    R.append(r.real + 1j*((r.imag + np.pi) % (2*np.pi) - np.pi))
    return np.array(R)

def newton(u, fixed, iters=80):
    n = len(u)
    for it in range(iters):
        F = residual(u, fixed); nF = np.linalg.norm(F)
        if not np.isfinite(nF): return None
        if nF < 1e-12: return u
        J = np.zeros((n, n), complex); h = 1e-7
        for k in range(n):
            du = np.zeros(n, complex); du[k] = h
            J[:, k] = (residual(u+du, fixed) - residual(u-du, fixed))/(2*h)
        try: du = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError: return None
        st = 1.0
        while st > 1e-9:
            F1 = residual(u+st*du, fixed)
            if np.all(np.isfinite(F1)) and np.linalg.norm(F1) < nF*(1-0.05*st): break
            st /= 2
        if st <= 1e-9: return None
        u = u + st*du
    return None

if __name__ == '__main__':
    nst = int(sys.argv[1])
    rng0 = np.random.default_rng(99)
    fixed = rng0.normal(size=3) + 1j*rng0.normal(size=3)
    ok = 0
    sols = []
    for seed in range(nst):
        rng = np.random.default_rng(seed)
        u = rng.normal(size=7) + 1j*rng.normal(size=7)
        r = newton(u, fixed)
        if r is not None:
            ok += 1; sols.append(r)
    print('converged', ok, 'of', nst)
    # distinct solutions (up to obvious symmetries) by z-values of t[:,1]
    zs = []
    for r in sols:
        key = np.round(np.sort_complex((r[0:3]**2 + r[0:3]**-2)), 6)
        if not any(np.allclose(key, k2) for k2 in zs): zs.append(key)
    print('distinct solutions (by z of second orbits):', len(zs))
    np.save('klein6t_sols.npy', np.array(sols)); np.save('klein6t_fixed.npy', fixed)
