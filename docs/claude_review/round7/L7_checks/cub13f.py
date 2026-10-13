"""NUMERICAL (complex float) twin of cub13q.build, for exploration only."""
import numpy as np
import itertools

P = np.polynomial.polynomial  # low-first coefficient arrays


def lin(r):
    return np.array([-r, 1.0], dtype=complex)


def build(tau, s34, s35, s45):
    A = 1 / ((tau - 1) * (tau - s45))
    B = 1 / (tau * (tau - 1) * (tau - s35))
    G = 1 / (tau * (tau - s34))
    c2 = A * P.polymul(lin(1), lin(s45))
    c3 = B * P.polymul(P.polymul(lin(0), lin(1)), lin(s35))
    c4 = G * P.polymul(lin(0), lin(s34))
    s36 = (B * s35 + G * s34) / (B * tau)
    K = P.polyval(s36, c3) / ((s36 - s34) * (s36 - s35) * (s36 - s45))
    c1 = K * P.polymul(P.polymul(lin(s34), lin(s35)), lin(s45))
    s246 = A * s45 / ((A - G) * tau)
    s256 = (A * s45 / B) / tau
    s346 = s35 * s45 / s36
    s356 = K * s34 * s45 / ((K - B) * s36)
    xs = np.array([s246, s256, s346, s356])
    ys = np.array([P.polyval(s246, c2), P.polyval(s256, c2), P.polyval(s346, c1), P.polyval(s356, c1)])
    V = np.vander(xs, 4, increasing=True)
    c0 = np.linalg.solve(V, ys)
    R, rem = P.polydiv(P.polysub(c1, c2), lin(s45))
    T = P.polysub(c0, c1)
    return dict(c=[c0, c1, c2, c3, c4], R=R, T=T, s36=s36, K=K, xs=xs)


def residual(v, tau, s34, s35):
    s45, sig = v
    d = build(tau, s34, s35, s45)
    return np.array([P.polyval(sig, d['R']), P.polyval(sig, d['T'])])


def newton(v, tau, s34, s35, it=60):
    v = np.array(v, dtype=complex)
    for _ in range(it):
        f = residual(v, tau, s34, s35)
        h = 1e-7
        J = np.zeros((2, 2), complex)
        for k in range(2):
            dv = np.zeros(2, complex)
            dv[k] = h
            J[:, k] = (residual(v + dv, tau, s34, s35) - residual(v - dv, tau, s34, s35)) / (2 * h)
        try:
            step = np.linalg.solve(J, -f)
        except np.linalg.LinAlgError:
            return None
        v = v + step
        if np.linalg.norm(step) < 1e-13 * (1 + np.linalg.norm(v)):
            break
    if np.linalg.norm(residual(v, tau, s34, s35)) > 1e-9:
        return None
    return v


# --- collision analysis in the y-chart: y_{k+1} = c_k, y_6 = 0, label 7 at infinity -----------
def yforms(c):
    """return the 6 y's as cubic forms (low-first length 4 arrays)"""
    ys = [np.concatenate([np.array(ck, complex), np.zeros(4 - len(ck))]) for ck in c]
    ys.append(np.zeros(4, complex))
    return ys  # index 0..5 -> labels 1..6


def events(c, tol=1e-6):
    """points s in P^1 (s = inf allowed) where some y's collide; returns list of (s, clusters)."""
    ys = yforms(c)
    pts = []
    for i, j in itertools.combinations(range(6), 2):
        d = ys[i] - ys[j]
        dd = np.trim_zeros(d, 'b')
        if len(dd) < 4:
            pts.append(np.inf)  # root at infinity (as a cubic form)
        if len(dd) > 1:
            pts += list(P.polyroots(dd))
    uniq = []
    for p in pts:
        if all(not ((p is np.inf and q is np.inf) or (p is not np.inf and q is not np.inf and abs(p - q) < tol)) for q in uniq):
            uniq.append(p)
    ev = []
    for p in uniq:
        if p is np.inf:
            vals = [y[3] for y in ys]
        else:
            vals = [P.polyval(p, y) for y in ys]
        groups = []
        for i in range(6):
            for g in groups:
                if abs(vals[i] - vals[g[0]]) < 1e-6 * (1 + abs(vals[i])):
                    g.append(i)
                    break
            else:
                groups.append([i])
        big = [tuple(x + 1 for x in g) for g in groups if len(g) >= 2]
        ev.append((p, big))
    return ev

TRIPLES13 = {(4, 5, 6), (3, 5, 6), (3, 4, 6), (3, 4, 5), (2, 5, 6), (2, 4, 6), (2, 3, 6), (2, 4, 5),
             (1, 3, 5), (1, 3, 4), (1, 2, 5), (1, 2, 4), (1, 2, 3)}
