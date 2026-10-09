"""NUMERICAL: toric-chart class-beta system (newton_log.py) plus a log-barrier equation that
excludes coincidences of special points and clusters larger than prescribed."""
import numpy as np, itertools, sys
from newton_beta import L, Aps, Ts, Aidx, PS, eqs, special_points

def xval(c, t, i, s):
    idx = [Aidx[A] for A in Aps if i in A]
    return c[i-1]*np.prod(s - t[idx])

def F_all(c, t, s, om):
    R = []
    for (T, i, j) in eqs:
        x = s[Ts.index(T)]
        Pi = PS[(i, j)]; Pj = PS[(j, i)]
        r = (c[i-1]/c[j-1])*np.prod((x - t[Pi])/(x - t[Pj]))
        R.append(np.log(r))
    # barrier
    pts = np.concatenate([t, s])
    sm = 0
    for a in range(len(pts)):
        for b in range(a+1, len(pts)):
            sm = sm + np.log(pts[a] - pts[b])
    for n, T in enumerate(Ts):
        i0 = min(T)
        xi = xval(c, t, i0, s[n])
        for k in L:
            if k not in T:
                sm = sm + np.log(xval(c, t, k, s[n]) - xi) - np.log(xi)
    for a in range(5):
        for b in range(a+1, 5):
            sm = sm + np.log(c[a] - c[b])
    r = sm + om
    R.append(r)
    R = np.array(R)
    R = R.real + 1j*((R.imag + np.pi) % (2*np.pi) - np.pi)
    return R

def run(seed, iters=150, nfix=9):
    rng = np.random.default_rng(seed)
    cp = lambda n: rng.normal(size=n) + 1j*rng.normal(size=n)
    c = cp(5); c[0] = 1.0
    t = cp(30); t[0] = 0.0; t[1] = 1.0
    fixed = [0, 1] + list(rng.choice(np.arange(2, 30), nfix, replace=False))
    free_t = [k for k in range(30) if k not in fixed]
    s = cp(16); om = np.array([0j])
    def pack(): return np.concatenate([c[1:], t[free_t], s, om])
    def unpack(v):
        c[1:] = v[:4]; t[free_t] = v[4:4+len(free_t)]; s[:] = v[4+len(free_t):4+len(free_t)+16]; om[0] = v[-1]
    v = pack(); n = len(v)
    def Fv(v):
        unpack(v); return F_all(c, t, s, om[0])
    for it in range(iters):
        F = Fv(v); nF = np.linalg.norm(F)
        if not np.isfinite(nF): return None
        if nF < 1e-11:
            unpack(v); return c.copy(), t.copy(), s.copy()
        J = np.zeros((len(F), n), complex); h = 1e-7
        for k in range(n):
            dv = np.zeros(n, complex); dv[k] = h
            J[:, k] = (Fv(v+dv) - Fv(v-dv))/(2*h)
        dv = np.linalg.lstsq(J, -F, rcond=None)[0]
        st = 1.0
        while st > 1e-9:
            F1 = Fv(v + st*dv)
            if np.all(np.isfinite(F1)) and np.linalg.norm(F1) < nF*(1-0.05*st): break
            st /= 2
        if st <= 1e-9: return None
        v = v + st*dv
    return None

if __name__ == '__main__':
    nst = int(sys.argv[1]); start = int(sys.argv[2])
    for seed in range(start, start+nst):
        r = run(seed)
        if r is None:
            print(seed, 'fail', flush=True); continue
        c, t, s = r
        pts, msg = special_points(c, t, s)
        np.save('nb_sol_%d.npy' % seed, np.concatenate([c, t, s]))
        print(seed, 'CONVERGED', msg, flush=True)
