"""NUMERICAL: one-involution ansatz (genus 0) for class-beta curves in Mbar_{0,7}.
tau: t -> -t moves only a.  b, c, d are functions of y = t^2 of degree 4; a has degree 8 in t.
Unknowns: for each anchor x in {0,1,inf}: t[x][B''] for the 7 nonempty B'' of {b,c,d} (a=x there too) and
t[x]['a'] (a=x alone): 24; mover points t_bcd, t_bc, t_bd, t_cd: 4; kappas ka,kb,kc,kd: 4; omega: 1.
Expected solution dimension 3 (2 + scaling t -> lambda t)."""
import numpy as np, itertools, sys
SUBS = [frozenset(s) for r in (1, 2, 3) for s in itertools.combinations('bcd', r)]  # 7
NU = 24 + 4 + 4 + 1

def unpack(u):
    T = u[:24].reshape(3, 8)  # columns: 7 SUBS then 'a'
    tm = u[24:28]; k = u[28:32]; om = u[32]
    return T, tm, k, om

def polys(u):
    T, tm, k, om = unpack(u)
    PA = [np.poly(T[x]) for x in range(3)]
    PM = {}
    for m in 'bcd':
        cols = [n for n, S in enumerate(SUBS) if m in S]
        PM[m] = [np.poly(T[x, cols]**2) for x in range(3)]
    return T, tm, k, om, PA, PM

def fval(P, kk, w):
    return kk*np.polyval(P[0], w)/np.polyval(P[2], w)

def residual(u):
    T, tm, k, om, PA, PM = polys(u)
    R = []
    deps = [(PA, k[0])] + [(PM[m], k[i+1]) for i, m in enumerate('bcd')]
    for P, kk in deps:
        dep = kk*P[0] - P[2] - (kk-1)*P[1]
        R += list(dep[1:]/ (1 + np.abs(P[2][1:])))
    A = lambda t: fval(PA, k[0], t)
    M = {m: (lambda t, m=m, kk=k[i+1]: fval(PM[m], kk, t**2)) for i, m in enumerate('bcd')}
    tbcd, tbc, tbd, tcd = tm
    g = A(tbcd); R += [M['b'](tbcd) - g, M['c'](tbcd) - g, M['d'](tbcd) - g]
    g = A(tbc); R += [M['b'](tbc) - g, M['c'](tbc) - g]
    g = A(tbd); R += [M['b'](tbd) - g, M['d'](tbd) - g]
    g = A(tcd); R += [M['c'](tcd) - g, M['d'](tcd) - g]
    # barrier: all 28 special y=t^2 values distinct (also t != 0, inf)
    ys = list((T**2).ravel()) + list(tm**2)
    sm = 0
    for i in range(len(ys)):
        sm = sm + np.log(ys[i])
        for j in range(i+1, len(ys)):
            sm = sm + np.log(ys[i] - ys[j])
    r = sm + om
    R.append(r.real + 1j*((r.imag + np.pi) % (2*np.pi) - np.pi))
    # fix scaling t -> lambda t: T[0][0] = 1
    R.append(T[0, 0] - 1)
    return np.array(R)

def gn(u, iters=100, tol=1e-11):
    n = len(u)
    for it in range(iters):
        F = residual(u); nF = np.linalg.norm(F)
        if not np.isfinite(nF): return None, it
        if nF < tol: return u, it
        J = np.zeros((len(F), n), complex); h = 1e-7
        for kk in range(n):
            du = np.zeros(n, complex); du[kk] = h
            J[:, kk] = (residual(u+du) - residual(u-du))/(2*h)
        du = np.linalg.lstsq(J, -F, rcond=None)[0]
        st = 1.0
        while st > 1e-9:
            F1 = residual(u + st*du)
            if np.all(np.isfinite(F1)) and np.linalg.norm(F1) < nF*(1-0.05*st): break
            st /= 2
        if st <= 1e-9: return None, it
        u = u + st*du
    return None, iters

if __name__ == '__main__':
    nst = int(sys.argv[1]); start = int(sys.argv[2])
    for seed in range(start, start+nst):
        rng = np.random.default_rng(seed)
        u = rng.normal(size=NU) + 1j*rng.normal(size=NU)
        u, it = gn(u)
        print(seed, 'OK' if u is not None else 'fail', it, flush=True)
        if u is not None: np.save('inv1_sol_%d.npy' % seed, u)
