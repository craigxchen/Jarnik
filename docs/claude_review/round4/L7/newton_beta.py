"""NUMERICAL search for rational curves of class beta in Mbar_{0,7}.

Toric chart: heavy points 0, infinity; light points x_1..x_5 with
  x_i(t) = c_i * prod_{A' ni i} (t - t_{A'})       (A' nonempty proper subset of [5])
The 30 roots t_{A'} are the contacts with D_{0 u A'}.  For T subset [5], |T|>=3,
s_T is a common point where all x_i (i in T) agree; it must be a root of the
reduced differences R_ij = c_i prod_{A' ni i, not j}(t-t_A') - c_j prod_{A' ni j, not i}(t-t_A').
Unknowns: c (5), t (30), s (16).  Normalize c_1=1, t_{A0}=0, t_{A1}=1, and fix 9
further t's at random values -> 39 unknowns, 39 equations.
"""
import numpy as np, itertools, sys

L = [1, 2, 3, 4, 5]
Aps = [frozenset(c) for r in range(1, 5) for c in itertools.combinations(L, r)]  # 30
Ts = [frozenset(c) for r in range(3, 6) for c in itertools.combinations(L, r)]   # 16
Aidx = {A: k for k, A in enumerate(Aps)}

def setP(i, j):
    return [Aidx[A] for A in Aps if i in A and j not in A]

PS = {(i, j): setP(i, j) for i in L for j in L if i != j}
eqs = []
for T in Ts:
    tl = sorted(T)
    i0 = tl[0]
    for j in tl[1:]:
        eqs.append((T, i0, j))
assert len(eqs) == 39

def residual_and_jac(c, t, s, free_c, free_t):
    """c: array len 5 (index i-1); t: len 30; s: dict T->value (array len 16)."""
    nE = len(eqs)
    F = np.zeros(nE, complex)
    # columns: free_c, free_t, s (16)
    nv = len(free_c) + len(free_t) + len(Ts)
    J = np.zeros((nE, nv), complex)
    cpos = {k: n for n, k in enumerate(free_c)}
    tpos = {k: len(free_c) + n for n, k in enumerate(free_t)}
    spos = {T: len(free_c) + len(free_t) + n for n, T in enumerate(Ts)}
    for e, (T, i, j) in enumerate(eqs):
        x = s[Ts.index(T)]
        Pi = PS[(i, j)]; Pj = PS[(j, i)]
        di = x - t[Pi]; dj = x - t[Pj]
        P = np.prod(di); Q = np.prod(dj)
        F[e] = c[i-1]*P - c[j-1]*Q
        # d/ds
        J[e, spos[T]] = c[i-1]*P*np.sum(1/di) - c[j-1]*Q*np.sum(1/dj)
        if (i-1) in cpos: J[e, cpos[i-1]] += P
        if (j-1) in cpos: J[e, cpos[j-1]] += -Q
        for n, k in enumerate(Pi):
            if k in tpos: J[e, tpos[k]] += -c[i-1]*P/di[n]
        for n, k in enumerate(Pj):
            if k in tpos: J[e, tpos[k]] += c[j-1]*Q/dj[n]
    return F, J

def run(seed, iters=200, verbose=False):
    rng = np.random.default_rng(seed)
    cplx = lambda n: rng.normal(size=n) + 1j*rng.normal(size=n)
    c = cplx(5); c[0] = 1.0
    t = cplx(30)
    t[0] = 0.0; t[1] = 1.0
    fixed_t = [0, 1] + list(rng.choice(np.arange(2, 30), 9, replace=False))
    free_t = [k for k in range(30) if k not in fixed_t]
    free_c = [1, 2, 3, 4]
    s = cplx(16)
    def pack():
        return np.concatenate([c[free_c], t[free_t], s])
    def unpack(v):
        c[free_c] = v[:4]; t[free_t] = v[4:4+len(free_t)]; s[:] = v[4+len(free_t):]
    v = pack()
    lam = 1.0
    for it in range(iters):
        unpack(v)
        F, J = residual_and_jac(c, t, s, free_c, free_t)
        nF = np.linalg.norm(F)/ (1+np.linalg.norm(v))
        if verbose: print(it, nF)
        if nF < 1e-13:
            return True, c.copy(), t.copy(), s.copy(), it
        try:
            dv = np.linalg.lstsq(J, -F, rcond=None)[0]
        except Exception:
            return False, None, None, None, it
        # damped step with backtracking on residual norm
        step = 1.0
        n0 = np.linalg.norm(F)
        while step > 1e-6:
            v1 = v + step*dv
            unpack(v1)
            F1, _ = residual_and_jac(c, t, s, free_c, free_t)
            if np.linalg.norm(F1) < n0*(1-0.1*step):
                break
            step /= 2
        v = v1
        if step <= 1e-6:
            unpack(v)
            return False, None, None, None, it
    unpack(v)
    return False, None, None, None, iters

def special_points(c, t, s):
    pts = {}
    for k, A in enumerate(Aps):
        pts[('0', A)] = t[k]
    for n, T in enumerate(Ts):
        pts[('d', T)] = s[n]
    # pair roots: roots of R_ij not among s_T
    for i, j in itertools.combinations(L, 2):
        Pi = PS[(i, j)]; Pj = PS[(j, i)]
        pi = np.poly(t[Pi]) * c[i-1]
        pj = np.poly(t[Pj]) * c[j-1]
        r = np.roots(pi - pj)
        known = [s[Ts.index(T)] for T in Ts if i in T and j in T]
        rem = list(r)
        for kv in known:
            d = [abs(x - kv) for x in rem]
            m = int(np.argmin(d))
            if d[m] > 1e-6:
                return None, 'known root missing'
            rem.pop(m)
        assert len(rem) == 1
        pts[('d', frozenset([i, j]))] = rem[0]
    return pts, 'ok'

if __name__ == '__main__':
    nst = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    found = 0
    for seed in range(nst):
        ok, c, t, s, it = run(seed)
        if not ok:
            continue
        pts, msg = special_points(c, t, s)
        if pts is None:
            print(seed, 'converged but', msg); continue
        vals = np.array(list(pts.values()))
        D = np.abs(vals[:, None] - vals[None, :]) + np.eye(len(vals))*1e9
        mind = D.min()
        cd = min(abs(c[a]-c[b]) for a in range(5) for b in range(a+1, 5))
        print(seed, 'converged it', it, 'min dist special pts %.3e' % mind, 'min |c_i-c_j| %.3e' % cd,
              'min|c| %.3e' % np.min(np.abs(c)), 'max|pt| %.3e' % np.max(np.abs(vals)))
        found += 1
    print('converged', found, 'of', nst)
