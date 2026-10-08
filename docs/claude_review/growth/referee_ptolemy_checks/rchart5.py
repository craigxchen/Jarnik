r"""NUMERICAL exploration (not a proof): balanced rational curves in Mbar_{0,6} (k=6, m=5)
in the 'R-chart':  p_0 = infinity,  R_i(x) = r_i/(x - a_i) + q_i(x),  q_1 = 0, deg q_i = 2.
The four-row cuts [5]\{i} are the simple poles a_i; every triple T={i,j,k} needs one point
rho_T with R_i = R_j = R_k; pair roots are the remaining zero of the degree-4 numerator of
R_i - R_j.  Balanced  <=>  all 25 special points distinct and leading coefficients of q_i distinct
and nonzero.  Converting back: X_i = E(x) R_i(x) with E = prod_l (x - a_l) gives polynomial rows
of degree 7 with X_i - X_j = (lead_i - lead_j) prod_{T contains i,j} (x - rho_T).
Usage: python3 rchart5.py trials [real]
"""
import sys, itertools
import numpy as np

m = 5
rows = list(range(1, 6))
triples = [frozenset(T) for T in itertools.combinations(rows, 3)]

def Rval(i, x, r, a, q):
    v = r[i] / (x - a[i])
    if i > 1:
        v = v + np.polyval(q[i], x)
    return v

def dR(i, x, r, a, q):
    v = -r[i] / (x - a[i]) ** 2
    if i > 1:
        v = v + np.polyval(np.polyder(q[i]), x)
    return v

def run(rng, real):
    rnd = (lambda n: rng.normal(size=n)) if real else (lambda n: (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2))
    dt = float if real else complex
    a = dict(zip(rows, rnd(5) * 2)); r = dict(zip(rows, rnd(5)))
    q = {i: rnd(3) for i in rows[1:]}
    rho = dict(zip(triples, rnd(10) * 2))
    fixedT = triples[:2] if True else []
    perm = list(range(10)); rng.shuffle(perm)
    fixedT = [triples[perm[0]], triples[perm[1]]]
    freeT = [T for T in triples if T not in fixedT]
    conds = []
    for T in triples:
        i, j, k = sorted(T)
        conds += [(T, i, j), (T, j, k)]
    def pack():
        return np.concatenate([np.concatenate([q[i] for i in rows[1:]]), np.array([rho[T] for T in freeT], dtype=dt)])
    def unpack(z):
        for n, i in enumerate(rows[1:]):
            q[i] = z[3 * n:3 * n + 3]
        for n, T in enumerate(freeT):
            rho[T] = z[12 + n]
    def F():
        return np.array([Rval(i, rho[T], r, a, q) - Rval(j, rho[T], r, a, q) for T, i, j in conds], dtype=dt)
    def J():
        Jm = np.zeros((20, 12 + 8), dtype=dt)
        for e, (T, i, j) in enumerate(conds):
            x = rho[T]
            pw = np.array([x * x, x, 1.0])
            if i > 1:
                Jm[e, 3 * (i - 2):3 * (i - 2) + 3] += pw
            if j > 1:
                Jm[e, 3 * (j - 2):3 * (j - 2) + 3] -= pw
            if T in freeT:
                Jm[e, 12 + freeT.index(T)] = dR(i, x, r, a, q) - dR(j, x, r, a, q)
        return Jm
    z = pack()
    for it in range(300):
        unpack(z); res = F(); nr = np.linalg.norm(res)
        if nr < 1e-14:
            break
        step = np.linalg.lstsq(J(), -res, rcond=None)[0]
        t = 1.0
        while t > 1e-8:
            zn = z + t * step; unpack(zn)
            if np.linalg.norm(F()) < nr:
                break
            t /= 2
        z = zn
    unpack(z)
    nr = np.linalg.norm(F())
    # pair roots: numerator of R_i - R_j
    special = {frozenset([l for l in rows if l != i]): a[i] for i in rows}
    special.update(rho)
    for i, j in itertools.combinations(rows, 2):
        qi = q[i] if i > 1 else np.zeros(3)
        qj = q[j] if j > 1 else np.zeros(3)
        num = np.polyadd(np.polymul(np.polysub(qi, qj), np.polymul([1, -a[i]], [1, -a[j]])),
                         np.polysub(r[i] * np.array([1, -a[j]]), r[j] * np.array([1, -a[i]])))
        roots = np.roots(num)
        known = [rho[T] for T in triples if i in T and j in T]
        rem = list(roots)
        for kv in known:
            d = [abs(x - kv) for x in rem]
            n = int(np.argmin(d))
            if d[n] > 1e-6:
                return nr, None, None, ('pairfail', i, j, d[n], kv, roots, [rho[T] for T in triples])
            rem.pop(n)
        special[frozenset((i, j))] = rem[0]
    vals = list(special.values())
    sep = min(abs(x - y) for x, y in itertools.combinations(vals, 2))
    leads = [0.0] + [q[i][0] for i in rows[1:]]
    lsep = min(abs(x - y) for x, y in itertools.combinations(leads, 2))
    return nr, sep, lsep, (a, r, q, special)

if __name__ == "__main__":
    trials = int(sys.argv[1]); real = len(sys.argv) > 2 and sys.argv[2] == "real"
    rng = np.random.default_rng(77 + (1 if real else 0))
    good = 0
    for tr in range(trials):
        out = run(rng, real)
        nr, sep, lsep = out[0], out[1], out[2]
        ok = nr < 1e-10 and sep is not None and sep > 1e-6 and lsep > 1e-6
        good += ok
        print(f"trial {tr}: residual {nr:.1e}, min sep {sep}, lead sep {lsep} -> {'NONDEGENERATE' if ok else 'fail'}")
    print("nondegenerate:", good, "/", trials)
