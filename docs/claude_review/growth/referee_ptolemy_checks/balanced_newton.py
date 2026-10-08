"""NUMERICAL exploration (not a proof): balanced degree-one rational curves in Mbar_{0,m+1}.

X-chart: anchor p_0 = infinity, rows X_1 = 0, X_2..X_m polynomials of degree D = 2^(m-2)-1
(the full block T=[m] is invisible in M_{0,m+1} and is dropped).  For every boundary pattern
T (subset of [m], 3 <= |T| <= m-1) there is one root rho_T with X_i(rho_T) = X_j(rho_T) for all
i,j in T.  For a pair {i,j}, X_i - X_j then has the D-1 roots rho_T (|T|>=3, T contains i,j);
its last root is the pair root rho_ij.  If all 2^m-m-2 roots are distinct and the leading
coefficients (0 for X_1) are distinct,
     X_i - X_j = (lead_i - lead_j) prod_{T contains i,j, T != [m]} (x - rho_T),
which is exactly a balanced curve (each boundary divisor met once, transversally).

Fixings: leading coefficients fixed to random distinct values; a random set of rho_T fixed.
Newton (least squares) in complex or real arithmetic.
Usage: python3 balanced_newton.py m trials [real]
"""
import sys, itertools
import numpy as np

def setup(m):
    rows = list(range(1, m + 1))
    big = [frozenset(T) for r in range(3, m) for T in itertools.combinations(rows, r)]
    D = 2 ** (m - 2) - 1
    conds = []
    for T in big:
        Ts = sorted(T)
        for a, b in zip(Ts, Ts[1:]):
            conds.append((T, a, b))
    return rows, big, D, conds

def X(coef, i, x):
    return 0.0 * x if i == 1 else np.polyval(coef[i], x)

def dX(coef, i, x):
    return 0.0 * x if i == 1 else np.polyval(np.polyder(coef[i]), x)

def pair_roots(m, rows, big, coef, rho):
    out = {}
    for i, j in itertools.combinations(rows, 2):
        li = coef[i][0] if i > 1 else 0.0
        lj = coef[j][0] if j > 1 else 0.0
        known = [rho[T] for T in big if i in T and j in T]
        x0 = 0.31 + 0.17j
        val = X(coef, i, x0) - X(coef, j, x0)
        out[frozenset((i, j))] = x0 - val / ((li - lj) * np.prod([x0 - r for r in known]))
    return out

def verify(m, rows, coef, allrho):
    worst = 0.0
    for i, j in itertools.combinations(rows, 2):
        li = coef[i][0] if i > 1 else 0.0
        lj = coef[j][0] if j > 1 else 0.0
        for x in (0.3 + 0.7j, -1.1 + 0.2j, 2.0 - 0.5j, 0.05 - 0.9j):
            lhs = X(coef, i, x) - X(coef, j, x)
            rhs = (li - lj) * np.prod([x - r for T, r in allrho.items() if i in T and j in T])
            worst = max(worst, abs(lhs - rhs) / (1 + abs(lhs)))
    return worst

def solve(m, rng, real=False, iters=200, scale=1.0):
    rows, big, D, conds = setup(m)
    nb = len(big); E = len(conds)
    ncoef = (m - 1) * D
    nfix = max(0, min(nb, ncoef + nb - E))
    rnd = (lambda n: rng.normal(size=n)) if real else (lambda n: (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2))
    leads = dict(zip(rows[1:], rnd(m - 1)))
    order = list(range(nb)); rng.shuffle(order)
    fixed = set(order[:nfix])
    rhov = rnd(nb) * scale
    rho = {T: rhov[n] for n, T in enumerate(big)}
    coef = {i: np.concatenate([[leads[i]], rnd(D)]) for i in rows[1:]}
    free = [T for n, T in enumerate(big) if n not in fixed]
    dt = float if real else complex
    def pack():
        return np.concatenate([np.concatenate([coef[i][1:] for i in rows[1:]]), np.array([rho[T] for T in free], dtype=dt)])
    def unpack(z):
        for n, i in enumerate(rows[1:]):
            coef[i] = np.concatenate([[leads[i]], z[n * D:(n + 1) * D]])
        for n, T in enumerate(free):
            rho[T] = z[ncoef + n]
    def F():
        return np.array([X(coef, a, rho[T]) - X(coef, b, rho[T]) for T, a, b in conds], dtype=dt)
    fpos = {T: n for n, T in enumerate(free)}
    def J():
        Jm = np.zeros((E, ncoef + len(free)), dtype=dt)
        for e, (T, a, b) in enumerate(conds):
            x = rho[T]
            powers = x ** np.arange(D - 1, -1, -1)
            for i, sgn in ((a, 1), (b, -1)):
                if i > 1:
                    n = i - 2
                    Jm[e, n * D:(n + 1) * D] += sgn * powers
            if T in fpos:
                Jm[e, ncoef + fpos[T]] = dX(coef, a, x) - dX(coef, b, x)
        return Jm
    z = pack()
    for it in range(iters):
        unpack(z)
        r = F(); nr = np.linalg.norm(r)
        if nr < 1e-13:
            break
        step = np.linalg.lstsq(J(), -r, rcond=None)[0]
        t = 1.0
        while t > 1e-8:
            zn = z + t * step
            unpack(zn)
            if np.linalg.norm(F()) < nr:
                break
            t /= 2
        z = zn
    unpack(z)
    nr = np.linalg.norm(F())
    pr = pair_roots(m, rows, big, coef, rho)
    allrho = dict(rho); allrho.update(pr)
    vals = list(allrho.values())
    sep = min(abs(a - b) for a, b in itertools.combinations(vals, 2))
    lsep = min(abs(x - y) for x, y in itertools.combinations([0.0] + list(leads.values()), 2))
    return nr, sep, lsep, coef, allrho, rows

if __name__ == "__main__":
    m = int(sys.argv[1]); trials = int(sys.argv[2]); real = len(sys.argv) > 3 and sys.argv[3] == "real"
    rng = np.random.default_rng(2026 + m + (100 if real else 0))
    rows, big, D, conds = setup(m)
    print(f"m={m} (k={m+1}): D={D}, blocks |T|>=3: {len(big)}, conditions {len(conds)}, "
          f"coef unknowns {(m-1)*D}, fixed roots {max(0, min(len(big), (m-1)*D+len(big)-len(conds)))}, "
          f"excess unknowns {(m-1)*D + len(big) - len(conds)}")
    good = 0
    for tr in range(trials):
        nr, sep, lsep, coef, allrho, rows = solve(m, rng, real=real)
        ok = nr < 1e-9 and sep > 1e-6 and lsep > 1e-6
        w = verify(m, rows, coef, allrho) if ok else float('nan')
        good += ok
        print(f"trial {tr}: residual {nr:.1e}, min sep of all {len(allrho)} roots {sep:.2e}, identity err {w:.1e} -> {'NONDEGENERATE' if ok else 'fail'}")
    print("nondegenerate solutions:", good, "of", trials)
