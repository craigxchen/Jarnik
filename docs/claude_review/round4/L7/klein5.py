"""NUMERICAL: Klein ansatz k=7, formulation with explicit numerator/denominator coefficients.
Unknowns: t[3][4], tcd, tabc, tabd, A=(NA,DA) deg 4 in v, B=(NB,DB) deg 4 in y, C=(NC,DC), D=(ND,DD) deg 2 in z,
plus barrier omega.  Each pair (N,D) is normalized by a fixed random linear form = 1."""
import numpy as np, sys

rng0 = np.random.default_rng(777)
NORM = {k: rng0.normal(size=n) + 1j*rng0.normal(size=n) for k, n in [('A', 10), ('B', 10), ('C', 6), ('D', 6)]}
sizes = [('t', 12), ('tcd', 1), ('tabc', 1), ('tabd', 1), ('A', 10), ('B', 10), ('C', 6), ('D', 6), ('om', 1)]
offs = {}; o = 0
for k, n in sizes: offs[k] = (o, o+n); o += n
NU = o

def get(u, k):
    a, b = offs[k]; return u[a:b]

def ND(u, k, deg):
    c = get(u, k); return c[:deg+1], c[deg+1:]

def ev(c, w):
    return np.polyval(c, w)

def residual(u, stage=3):
    t = get(u, 't').reshape(3, 4)
    v = t + 1/t; y = t**2; z = y + 1/y
    NA, DA = ND(u, 'A', 4); NB, DB = ND(u, 'B', 4); NC, DC = ND(u, 'C', 2); NDd, DDd = ND(u, 'D', 2)
    R = []
    for k in 'ABCD':
        R.append(np.dot(NORM[k], get(u, k)) - 1)
    def fib(N, Dn, w, x):
        if x == 0: return ev(N, w)
        if x == 1: return ev(N, w) - ev(Dn, w)
        return ev(Dn, w)
    sc = lambda w: 1.0/(1+abs(w))**2
    for x in range(3):
        for i in range(4):
            R.append(fib(NA, DA, v[x, i], x)*sc(v[x, i])**2)
            R.append(fib(NB, DB, y[x, i], x)*sc(y[x, i])**2)
        for i in (0, 1):
            R.append(fib(NC, DC, z[x, i], x)*sc(z[x, i]))
        for i in (1, 2):
            R.append(fib(NDd, DDd, z[x, i], x)*sc(z[x, i]))
    zs = list(z.ravel())
    def cross(N1, D1, w1, N2, D2, w2):
        return ev(N1, w1)*ev(D2, w2) - ev(N2, w2)*ev(D1, w1)
    if stage >= 1:
        s = get(u, 'tcd')[0]; vv, yy, zz = s + 1/s, s**2, s**2 + s**-2
        R += [cross(NC, DC, zz, NDd, DDd, zz), cross(NA, DA, vv, NC, DC, zz), cross(NB, DB, yy, NC, DC, zz)]
        zs.append(zz)
    if stage >= 2:
        s = get(u, 'tabc')[0]; vv, yy, zz = s + 1/s, s**2, s**2 + s**-2
        R += [cross(NA, DA, vv, NC, DC, zz), cross(NB, DB, yy, NC, DC, zz)]
        zs.append(zz)
    if stage >= 3:
        s = get(u, 'tabd')[0]; vv, yy, zz = s + 1/s, s**2, s**2 + s**-2
        R += [cross(NA, DA, vv, NDd, DDd, zz), cross(NB, DB, yy, NDd, DDd, zz)]
        zs.append(zz)
    sm = 0
    for i in range(len(zs)):
        sm = sm + np.log(zs[i] - 2) + np.log(zs[i] + 2)
        for j in range(i+1, len(zs)):
            sm = sm + np.log(zs[i] - zs[j])
    r = sm + get(u, 'om')[0]
    R.append(r.real + 1j*((r.imag + np.pi) % (2*np.pi) - np.pi))
    return np.array(R)

def gn(u, stage, iters=80, tol=1e-12):
    n = len(u)
    for it in range(iters):
        F = residual(u, stage); nF = np.linalg.norm(F)
        if not np.isfinite(nF): return None
        if nF < tol: return u
        J = np.zeros((len(F), n), complex); h = 1e-7
        for k in range(n):
            du = np.zeros(n, complex); du[k] = h
            J[:, k] = (residual(u+du, stage) - residual(u-du, stage))/(2*h)
        du = np.linalg.lstsq(J, -F, rcond=None)[0]
        st = 1.0
        while st > 1e-9:
            u1 = u + st*du; F1 = residual(u1, stage)
            if np.all(np.isfinite(F1)) and np.linalg.norm(F1) < nF*(1-0.05*st): break
            st /= 2
        if st <= 1e-9: return None
        u = u1
    return None

if __name__ == '__main__':
    nst = int(sys.argv[1]); start = int(sys.argv[2])
    for seed in range(start, start+nst):
        rng = np.random.default_rng(seed)
        u = rng.normal(size=NU) + 1j*rng.normal(size=NU)
        ok = True
        for stage in range(4):
            u = gn(u, stage)
            if u is None:
                print(seed, 'failed at stage', stage, flush=True); ok = False; break
        if ok:
            np.save('k5sol_%d.npy' % seed, u); print(seed, 'SOLVED', flush=True)
