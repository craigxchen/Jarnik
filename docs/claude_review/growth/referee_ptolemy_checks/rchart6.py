r"""NUMERICAL exploration (not a proof): balanced rational curves in Mbar_{0,7} (k=7, m=6).

R-chart: p_0 = infinity, R_i(x) = r_i/(x-a_i) + q_i(x), deg q_i <= 9, q_1 = 0, r_1 = 1.
The six 5-row cuts are the poles a_i.  Every T with |T| in {3,4} needs a point rho_T with
R_i(rho_T) all equal for i in T; pair points are the remaining zeros.  Balanced iff all 56
special points are distinct, all r_i != 0 and the leading coefficients of q_i (0 for q_1) are
distinct.  A set of triple points is fixed (gauge + parameters); the system is square.
Usage: python3 rchart6.py trials seed
"""
import sys, itertools
import numpy as np

rows = list(range(1, 7))
triples = [frozenset(T) for T in itertools.combinations(rows, 3)]
quads = [frozenset(T) for T in itertools.combinations(rows, 4)]
DEG = 9

def choose_fixed(rng, nfix=11):
    # prefer a set of triples meeting every 4-set at least twice
    best = None
    for _ in range(2000):
        F = set(rng.choice(len(triples), nfix, replace=False).tolist())
        cov = min(sum(1 for n in F if triples[n] <= Q) for Q in quads)
        if best is None or cov > best[0]:
            best = (cov, F)
        if cov >= 2:
            break
    return [triples[n] for n in sorted(best[1])]

def build(rng):
    c = lambda n: (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
    fixedT = choose_fixed(rng)
    fixedv = dict(zip(fixedT, c(len(fixedT)) * 2))
    freeB = [T for T in triples if T not in fixedT] + quads
    conds = []
    for T in triples + quads:
        Ts = sorted(T)
        for a, b in zip(Ts, Ts[1:]):
            conds.append((T, a, b))
    # unknown layout: a_1..a_6 (6), r_2..r_6 (5), q_2..q_6 (5*(DEG+1)), free block points
    nq = DEG + 1
    idx_a = {i: i - 1 for i in rows}
    idx_r = {i: 6 + (i - 2) for i in rows[1:]}
    idx_q = {i: 11 + (i - 2) * nq for i in rows[1:]}
    idx_p = {T: 11 + 5 * nq + n for n, T in enumerate(freeB)}
    nz = 11 + 5 * nq + len(freeB)
    assert nz == len(conds), (nz, len(conds))
    z0 = np.concatenate([c(6) * 2, c(5), c(5 * nq) * 0.3, c(len(freeB)) * 2])
    return dict(fixedv=fixedv, freeB=freeB, conds=conds, idx_a=idx_a, idx_r=idx_r, idx_q=idx_q,
                idx_p=idx_p, nz=nz, z0=z0, nq=nq)

def point(S, z, T):
    return z[S['idx_p'][T]] if T in S['idx_p'] else S['fixedv'][T]

def Rv(S, z, i, x):
    a = z[S['idx_a'][i]]
    r = 1.0 if i == 1 else z[S['idx_r'][i]]
    v = r / (x - a)
    if i > 1:
        st = S['idx_q'][i]
        v = v + np.polyval(z[st:st + S['nq']], x)
    return v

def dRv(S, z, i, x):
    a = z[S['idx_a'][i]]
    r = 1.0 if i == 1 else z[S['idx_r'][i]]
    v = -r / (x - a) ** 2
    if i > 1:
        st = S['idx_q'][i]
        v = v + np.polyval(np.polyder(z[st:st + S['nq']]), x)
    return v

def F(S, z):
    return np.array([Rv(S, z, a, point(S, z, T)) - Rv(S, z, b, point(S, z, T)) for T, a, b in S['conds']])

def J(S, z):
    Jm = np.zeros((len(S['conds']), S['nz']), dtype=complex)
    nq = S['nq']
    for e, (T, a, b) in enumerate(S['conds']):
        x = point(S, z, T)
        for i, s in ((a, 1), (b, -1)):
            ai = z[S['idx_a'][i]]
            r = 1.0 if i == 1 else z[S['idx_r'][i]]
            Jm[e, S['idx_a'][i]] += s * r / (x - ai) ** 2
            if i > 1:
                Jm[e, S['idx_r'][i]] += s / (x - ai)
                st = S['idx_q'][i]
                Jm[e, st:st + nq] += s * x ** np.arange(nq - 1, -1, -1)
        if T in S['idx_p']:
            Jm[e, S['idx_p'][T]] = dRv(S, z, a, x) - dRv(S, z, b, x)
    return Jm

def newton(S, z, iters=300):
    for it in range(iters):
        f = F(S, z); nf = np.linalg.norm(f)
        if nf < 1e-13:
            break
        step = np.linalg.lstsq(J(S, z), -f, rcond=None)[0]
        t = 1.0
        while t > 1e-9:
            zn = z + t * step
            if np.linalg.norm(F(S, zn)) < nf:
                break
            t /= 2
        z = zn
    return z, np.linalg.norm(F(S, z))

def analyse(S, z):
    """Return ((sep, lead_sep, rmin, special, sep_big, rem), None) or (None, reason).  Pair points
    are obtained by polynomial division of the numerator of R_i - R_j by the known factors."""
    special = {}
    for i in rows:
        special[frozenset(r for r in rows if r != i)] = z[S['idx_a'][i]]
    for T in triples + quads:
        special[T] = point(S, z, T)
    big = list(special.values())
    sep_big = min(abs(x - y) for x, y in itertools.combinations(big, 2))
    worst_rem = 0.0
    for i, j in itertools.combinations(rows, 2):
        ai, aj = z[S['idx_a'][i]], z[S['idx_a'][j]]
        ri = 1.0 if i == 1 else z[S['idx_r'][i]]
        rj = 1.0 if j == 1 else z[S['idx_r'][j]]
        qi = np.zeros(S['nq']) if i == 1 else z[S['idx_q'][i]:S['idx_q'][i] + S['nq']]
        qj = np.zeros(S['nq']) if j == 1 else z[S['idx_q'][j]:S['idx_q'][j] + S['nq']]
        num = np.polyadd(np.polymul(np.polysub(qi, qj), np.polymul([1, -ai], [1, -aj])),
                         np.polysub(ri * np.array([1, -aj]), rj * np.array([1, -ai])))
        known = [special[T] for T in triples + quads if i in T and j in T]
        den = np.poly(known)
        quo, rem = np.polydiv(num, den)
        worst_rem = max(worst_rem, np.max(np.abs(rem)) / (1 + np.max(np.abs(num))))
        if len(quo) != 2 or abs(quo[0]) < 1e-9:
            return None, f'pair {i}{j}: degree drop'
        special[frozenset((i, j))] = -quo[1] / quo[0]
    vals = list(special.values())
    sep = min(abs(x - y) for x, y in itertools.combinations(vals, 2))
    leads = [0.0] + [z[S['idx_q'][i]] for i in rows[1:]]
    lsep = min(abs(x - y) for x, y in itertools.combinations(leads, 2))
    rmin = min([1.0] + [abs(z[S['idx_r'][i]]) for i in rows[1:]])
    return (sep, lsep, rmin, special, sep_big, worst_rem), None

if __name__ == "__main__":
    trials = int(sys.argv[1]); seed = int(sys.argv[2])
    rng = np.random.default_rng(seed)
    good = 0
    for tr in range(trials):
        S = build(rng)
        z, nf = newton(S, S['z0'])
        if nf > 1e-9:
            print(f"trial {tr}: no convergence ({nf:.1e})", flush=True); continue
        out, why = analyse(S, z)
        np.save(f"r6sol_seed{seed}_tr{tr}.npy", z)
        if out is None:
            print(f"trial {tr}: converged; {why}", flush=True); continue
        sep, lsep, rmin, special, sep_big, rem = out
        ok = sep > 1e-5 and lsep > 1e-5 and rmin > 1e-5 and rem < 1e-8
        good += ok
        print(f"trial {tr}: converged; min sep {sep:.2e} (big {sep_big:.2e}), lead sep {lsep:.2e}, min|r| {rmin:.2e}, rem {rem:.1e} -> {'NONDEGENERATE' if ok else 'degenerate'}", flush=True)
    print("nondegenerate:", good, "/", trials)
