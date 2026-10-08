r"""NUMERICAL exploration (not a proof), k=6 (m=5), R-chart with ALL ten triple points fixed:
R_i(x) = r_i/(x - a_i) + q_i(x), deg q_i <= 2, q_1 = 0 (common quadratic shift gauge), r_1 = 1
(scaling gauge), a_1 fixed.  Unknowns: a_2..a_5, r_2..r_5, q_2..q_5 (20); equations: for every
triple T={i<j<k} with fixed point rho_T: R_i = R_j and R_j = R_k at rho_T (20).
A solution is balanced iff the 5 poles, 10 triple points and 10 pair points are distinct,
all r_i != 0, and the leading coefficients 0, lead(q_2..q_5) are distinct.
"""
import sys, itertools
import numpy as np

rows = [1, 2, 3, 4, 5]
triples = [frozenset(T) for T in itertools.combinations(rows, 3)]
conds = []
for T in triples:
    i, j, k = sorted(T)
    conds += [(T, i, j), (T, j, k)]

def unpack(z, a1):
    a = {1: a1}; r = {1: 1.0}; q = {1: np.zeros(3)}
    for n, i in enumerate(rows[1:]):
        a[i] = z[n]; r[i] = z[4 + n]; q[i] = z[8 + 3 * n: 11 + 3 * n]
    return a, r, q

def Rv(i, x, a, r, q):
    return r[i] / (x - a[i]) + np.polyval(q[i], x)

def F(z, a1, rho):
    a, r, q = unpack(z, a1)
    return np.array([Rv(i, rho[T], a, r, q) - Rv(j, rho[T], a, r, q) for T, i, j in conds])

def Jac(z, a1, rho):
    a, r, q = unpack(z, a1)
    J = np.zeros((20, 20), dtype=complex)
    for e, (T, i, j) in enumerate(conds):
        x = rho[T]
        for l, s in ((i, 1), (j, -1)):
            if l == 1:
                continue
            n = l - 2
            J[e, n] += s * r[l] / (x - a[l]) ** 2          # d/da_l
            J[e, 4 + n] += s / (x - a[l])                   # d/dr_l
            J[e, 8 + 3 * n: 11 + 3 * n] += s * np.array([x * x, x, 1.0])
    return J

def newton(z, a1, rho, iters=200):
    for it in range(iters):
        f = F(z, a1, rho); nf = np.linalg.norm(f)
        if nf < 1e-14:
            break
        step = np.linalg.lstsq(Jac(z, a1, rho), -f, rcond=None)[0]
        t = 1.0
        while t > 1e-9:
            zn = z + t * step
            if np.linalg.norm(F(zn, a1, rho)) < nf:
                break
            t /= 2
        z = zn
    return z, np.linalg.norm(F(z, a1, rho))

def analyse(z, a1, rho):
    a, r, q = unpack(z, a1)
    special = {('pole', i): a[i] for i in rows}
    special.update({('tri', T): v for T, v in rho.items()})
    for i, j in itertools.combinations(rows, 2):
        num = np.polyadd(np.polymul(np.polysub(q[i], q[j]), np.polymul([1, -a[i]], [1, -a[j]])),
                         np.polysub(r[i] * np.array([1, -a[j]]), r[j] * np.array([1, -a[i]])))
        num = np.trim_zeros(num, 'f')
        roots = list(np.roots(num))
        for T in triples:
            if i in T and j in T:
                d = [abs(x - rho[T]) for x in roots]
                n = int(np.argmin(d))
                if d[n] > 1e-6:
                    return None, 'pair numerator misses a triple point'
                roots.pop(n)
        if len(roots) != 1:
            return None, f'pair {i}{j}: numerator degree problem ({len(roots)} leftover)'
        special[('pair', (i, j))] = roots[0]
    vals = list(special.values())
    sep = min(abs(x - y) for x, y in itertools.combinations(vals, 2))
    leads = [0.0] + [q[i][0] for i in rows[1:]]
    lsep = min(abs(x - y) for x, y in itertools.combinations(leads, 2))
    rmin = min(abs(r[i]) for i in rows)
    return (sep, lsep, rmin), special

if __name__ == "__main__":
    trials = int(sys.argv[1])
    rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    good = []
    for tr in range(trials):
        c = lambda n: (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
        rho = dict(zip(triples, c(10) * 2))
        a1 = c(1)[0] * 2
        z0 = np.concatenate([c(4) * 2, c(4), c(12)])
        z, nf = newton(z0, a1, rho)
        if nf > 1e-10:
            print(f"trial {tr}: no convergence {nf:.1e}"); continue
        info, special = analyse(z, a1, rho)
        if info is None:
            print(f"trial {tr}: converged, {special}"); continue
        sep, lsep, rmin = info
        ok = sep > 1e-6 and lsep > 1e-6 and rmin > 1e-6
        print(f"trial {tr}: converged; min sep {sep:.2e}, lead sep {lsep:.2e}, min|r| {rmin:.2e} -> {'NONDEGENERATE' if ok else 'degenerate'}")
        if ok:
            good.append((z, a1, rho))
    print("nondegenerate:", len(good))
