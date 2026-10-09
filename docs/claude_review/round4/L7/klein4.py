"""NUMERICAL: staged Gauss-Newton for the Klein ansatz (k=7), t-coordinates as in klein3.py."""
import numpy as np, sys
from klein3 import unpack, funcs

def residual_stage(u, stage):
    t, tcd, tabc, tabd, kA, kB, kC, kD, om = unpack(u)
    t_, v, y, z, PA, PB, QC, QD, A, B, C, D = funcs(u)
    R = []
    for P, k in [(PA, kA), (PB, kB), (QC, kC), (QD, kD)]:
        dep = k*P[0] - P[2] - (k-1)*P[1]
        R += list(dep[1:])
    vyz = lambda s: (s + 1/s, s**2, s**2 + s**-2)
    zs = list(z.ravel())
    if stage >= 1:
        vv, yy, zz = vyz(tcd); g = C(zz)
        R += [D(zz) - g, A(vv) - g, B(yy) - g]; zs.append(zz)
    if stage >= 2:
        vv, yy, zz = vyz(tabc); g = C(zz)
        R += [A(vv) - g, B(yy) - g]; zs.append(zz)
    if stage >= 3:
        vv, yy, zz = vyz(tabd); g = D(zz)
        R += [A(vv) - g, B(yy) - g]; zs.append(zz)
    s = 0
    for i in range(len(zs)):
        s = s + np.log(zs[i] - 2) + np.log(zs[i] + 2)
        for j in range(i+1, len(zs)):
            s = s + np.log(zs[i] - zs[j])
    r = s + om
    R.append(r.real + 1j*((r.imag + np.pi) % (2*np.pi) - np.pi))
    return np.array(R)

def gn(u, stage, iters=60, tol=1e-11):
    n = len(u)
    for it in range(iters):
        F = residual_stage(u, stage)
        nF = np.linalg.norm(F)
        if not np.isfinite(nF): return None
        if nF < tol: return u
        J = np.zeros((len(F), n), complex); h = 1e-7
        for k in range(n):
            du = np.zeros(n, complex); du[k] = h
            J[:, k] = (residual_stage(u+du, stage) - residual_stage(u-du, stage))/(2*h)
        du = np.linalg.lstsq(J, -F, rcond=None)[0]
        st = 1.0
        while st > 1e-9:
            u1 = u + st*du; F1 = residual_stage(u1, stage)
            if np.all(np.isfinite(F1)) and np.linalg.norm(F1) < nF*(1-0.05*st): break
            st /= 2
        if st <= 1e-9: return None
        u = u1
    return None

if __name__ == '__main__':
    nst = int(sys.argv[1]); start = int(sys.argv[2])
    for seed in range(start, start+nst):
        rng = np.random.default_rng(seed)
        u = rng.normal(size=20) + 1j*rng.normal(size=20)
        ok = True
        for stage in range(4):
            u = gn(u, stage)
            if u is None:
                print(seed, 'failed at stage', stage, flush=True); ok = False; break
        if ok:
            np.save('k4sol_%d.npy' % seed, u)
            print(seed, 'SOLVED', flush=True)
