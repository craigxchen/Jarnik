"""NUMERICAL geometric genus of an irreducible plane curve F(x,y)=0 via monodromy + Riemann-Hurwitz.
F given as sympy polynomial in symbols x, y."""
import sympy as sp, numpy as np, mpmath as mp
mp.mp.dps = 30

def monodromy_genus(F, x, y, verbose=False):
    P = sp.Poly(F, y)
    n = P.degree()
    coeffs = [sp.Poly(c, x) for c in P.all_coeffs()]  # leading first
    def roots_at(xv):
        cs = [complex(c.eval(xv)) if c.degree() >= 0 else 0j for c in coeffs]
        return np.roots(cs)
    disc = sp.Poly(sp.discriminant(P.as_expr(), y), x)
    lc = coeffs[0]
    crit = []
    for poly in (disc, lc):
        if poly.degree() > 0:
            for fac, mlt in sp.factor_list(poly.as_expr())[1]:
                fp = sp.Poly(fac, x)
                if fp.degree() <= 0: continue
                rts = fp.nroots(n=25, maxsteps=500)
                crit += [complex(r) for r in rts]
    # merge numerically equal
    pts = []
    for c in crit:
        if all(abs(c - p) > 1e-8*(1+abs(c)) for p in pts): pts.append(c)
    total = 0
    def track(center, rad, nsteps=4000):
        th = np.linspace(0, 2*np.pi, nsteps+1)
        cur = roots_at(center + rad)
        start = cur.copy()
        for k in range(1, nsteps+1):
            nxt = roots_at(center + rad*np.exp(1j*th[k]))
            # match
            perm = []
            used = set()
            for r in cur:
                dists = [abs(r - s) if j not in used else np.inf for j, s in enumerate(nxt)]
                j = int(np.argmin(dists)); used.add(j); perm.append(nxt[j])
            cur = np.array(perm)
        # permutation start -> cur
        pm = []
        for r in cur:
            pm.append(int(np.argmin([abs(r - s) for s in start])))
        return pm
    def ncycles(pm):
        seen = set(); c = 0
        for i in range(len(pm)):
            if i in seen: continue
            c += 1; j = i
            while j not in seen:
                seen.add(j); j = pm[j]
        return c
    for p in pts:
        others = [abs(p - q) for q in pts if q is not p]
        rad = min(others)/3 if others else 1.0
        pm = track(p, rad)
        if sorted(pm) != list(range(n)):
            pm = track(p, rad, 20000)
        cyc = ncycles(pm)
        total += n - cyc
        if verbose and cyc != n: print('branch', np.round(p, 6), 'cycles', cyc)
    R = 3*max(abs(p) for p in pts) + 3
    pm = track(0, R, 20000)
    cyc = ncycles(pm); total += n - cyc
    if verbose: print('infinity cycles', cyc)
    g = 1 - n + total/2
    return g, n, len(pts)
