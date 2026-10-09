"""NUMERICAL: log-ratio formulation of the class-beta system in the toric chart (see newton_beta.py)."""
import numpy as np, itertools, sys
from newton_beta import L, Aps, Ts, Aidx, PS, eqs, special_points

def FJ(c, t, s, free_c, free_t):
    nE = len(eqs); nv = len(free_c)+len(free_t)+len(Ts)
    F = np.zeros(nE, complex); J = np.zeros((nE, nv), complex)
    cpos = {k: n for n, k in enumerate(free_c)}
    tpos = {k: len(free_c)+n for n, k in enumerate(free_t)}
    off = len(free_c)+len(free_t)
    for e, (T, i, j) in enumerate(eqs):
        n = Ts.index(T); x = s[n]
        Pi = PS[(i, j)]; Pj = PS[(j, i)]
        di = x - t[Pi]; dj = x - t[Pj]
        r = (c[i-1]/c[j-1]) * np.prod(di/dj)
        F[e] = np.log(r)
        J[e, off+n] = np.sum(1/di) - np.sum(1/dj)
        if (i-1) in cpos: J[e, cpos[i-1]] += 1/c[i-1]
        if (j-1) in cpos: J[e, cpos[j-1]] += -1/c[j-1]
        for m, k in enumerate(Pi):
            if k in tpos: J[e, tpos[k]] += -1/di[m]
        for m, k in enumerate(Pj):
            if k in tpos: J[e, tpos[k]] += 1/dj[m]
    return F, J

def run(seed, iters=300, scale=1.0, nfix=9):
    rng = np.random.default_rng(seed)
    cplx = lambda n: (rng.normal(size=n) + 1j*rng.normal(size=n))*scale
    c = cplx(5); c[0] = 1.0
    t = cplx(30); t[0] = 0.0; t[1] = 1.0
    fixed_t = [0, 1] + list(rng.choice(np.arange(2, 30), nfix, replace=False))
    free_t = [k for k in range(30) if k not in fixed_t]
    free_c = [1, 2, 3, 4]
    s = cplx(16)
    def pack(): return np.concatenate([c[free_c], t[free_t], s])
    def unpack(v):
        c[free_c] = v[:4]; t[free_t] = v[4:4+len(free_t)]; s[:] = v[4+len(free_t):]
    v = pack()
    for it in range(iters):
        unpack(v)
        F, J = FJ(c, t, s, free_c, free_t)
        if not np.all(np.isfinite(F)) or not np.all(np.isfinite(J)):
            return False, None, None, None, it
        nF = np.linalg.norm(F)
        if nF < 1e-12:
            return True, c.copy(), t.copy(), s.copy(), it
        dv = np.linalg.lstsq(J, -F, rcond=None)[0]
        step = 1.0
        while step > 1e-8:
            v1 = v + step*dv
            unpack(v1)
            F1, _ = FJ(c, t, s, free_c, free_t)
            if np.all(np.isfinite(F1)) and np.linalg.norm(F1) < nF*(1-0.05*step):
                break
            step /= 2
        if step <= 1e-8:
            return False, None, None, None, it
        v = v1
    unpack(v)
    return False, None, None, None, iters

if __name__ == '__main__':
    nst = int(sys.argv[1]); start = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    good = []
    for seed in range(start, start+nst):
        ok, c, t, s, it = run(seed)
        if not ok: continue
        pts, msg = special_points(c, t, s)
        if pts is None:
            print(seed, 'conv but', msg); continue
        vals = np.array(list(pts.values()))
        D = np.abs(vals[:, None]-vals[None, :]) + np.eye(len(vals))*1e9
        cd = min(abs(c[a]-c[b]) for a in range(5) for b in range(a+1, 5))
        print(seed, 'conv it', it, 'minsep %.2e' % D.min(), 'min|ci-cj| %.2e' % cd, 'max|pt| %.2e' % np.abs(vals).max(), flush=True)
        good.append(seed)
    print('done', good)
