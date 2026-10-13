"""Separation lemma (SL) test on small slab-box models.  Exact integer arithmetic.

Model: d = n-1 coordinates k in [-T,T]^d; profile sign columns beta_q in {+-1}^n;
y_q(k) = <k, beta_q>/2 on the zero-sum slice = sum_{x<d} k_x (beta_{q,x}-beta_{q,n})/2  (integer).
Narrow slabs |y_q| <= s for q != p.  The wide prime p has level l(k) = y_p(k) (e_p = 1).
For k* in A, S = {k*} cup {k in A : l(k) > l(k*)}.
char0 separator: g in Q[k], deg < m, g = 0 on S minus k*, g(k*) = 1  <=>  rank Ev_S = rank Ev_{S\k*} + 1.
p-adic separator: g with g(k*) = 1, v_p(g(k)) >= l(k) - l(k*) on S \ k*.
No p-adic separator <=> exists integer c in V(S,m) with p-adically unique minimiser of v_p(c_k)+l(k) at k*.
"""
import itertools, random, sys
from math import comb
from fractions import Fraction

def monos(d, m):
    out = []
    for deg in range(m):
        for a in itertools.product(range(deg + 1), repeat=d):
            if sum(a) == deg:
                out.append(a)
    return out

def binom_val(k, a, T):
    v = 1
    for kx, ax in zip(k, a):
        v *= comb(kx + T, ax)
    return v

def rank_mod(rows, P):
    rows = [[x % P for x in r] for r in rows]
    rk = 0
    ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = None
        for i in range(rk, len(rows)):
            if rows[i][c]:
                piv = i; break
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], P - 2, P)
        rows[rk] = [(x * inv) % P for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c]:
                f = rows[i][c]
                rows[i] = [(x - f * y) % P for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def exact_rank(rows):
    # fraction-free (Bareiss-like) rank over Q using python ints
    M = [list(r) for r in rows]
    rk = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = None
        for i in range(rk, len(M)):
            if M[i][c] != 0:
                piv = i; break
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        pv = M[rk][c]
        for i in range(rk + 1, len(M)):
            if M[i][c]:
                f = M[i][c]
                M[i] = [pv * x - f * y for x, y in zip(M[i], M[rk])]
                # reduce by gcd to keep sizes small
                from math import gcd
                g = 0
                for x in M[i]:
                    g = gcd(g, x)
                if g > 1:
                    M[i] = [x // g for x in M[i]]
        rk += 1
    return rk

def padic_separator(S_pts, kstar_idx, levels, lstar, mon, T, p):
    """Return True if a p-adic separator exists.
    Columns: points of S. Rows: basis polys (binomial basis, integer valued).
    v_alpha = (g_alpha(k) * p^{-(l(k)-lstar)})_k.  L = Q-span cap Z_(p)^S. Separator exists iff some vector of
    (L mod p) is nonzero at k*.  We compute a Z_(p)-basis of L by p-saturation of the integer matrix
    B = (g_alpha(k) * p^{Lmax-(l(k)-lstar)}) (scaled by p^{Lmax}), then divide out.
    """
    Lmax = max(levels[i] - lstar for i in range(len(S_pts)))
    B = []
    for a in mon:
        B.append([binom_val(k, a, T) * p ** (Lmax - (levels[i] - lstar)) for i, k in enumerate(S_pts)])
    # p-saturation: we want the lattice L' = (Q-span of B) cap Z_(p)^S.  Algorithm: maintain integer rows
    # generating a full-rank (over Q) set; repeatedly find F_p-dependencies among rows mod p, replace a row
    # by (combination)/p.  Terminates when rows are independent mod p.  Work with Z_(p): any unit scalars ok.
    # First reduce to a Q-basis.
    rows = [r[:] for r in B]
    # Q-row-reduce to independent rows (use exact fraction-free elimination keeping integer rows)
    basis = []
    for r in rows:
        basis.append(r)
        if exact_rank(basis) < len(basis):
            basis.pop()
    rows = basis
    while True:
        # find dependency mod p
        n_r = len(rows)
        A = [[x % p for x in r] for r in rows]
        # augmented elimination tracking combinations
        comb_ = [[1 if i == j else 0 for j in range(n_r)] for i in range(n_r)]
        rk = 0
        ncol = len(A[0])
        piv_rows = []
        Awork = [a[:] for a in A]
        dep = None
        order = list(range(n_r))
        for c in range(ncol):
            piv = None
            for i in range(rk, n_r):
                if Awork[i][c]:
                    piv = i; break
            if piv is None:
                continue
            Awork[rk], Awork[piv] = Awork[piv], Awork[rk]
            comb_[rk], comb_[piv] = comb_[piv], comb_[rk]
            inv = pow(Awork[rk][c], p - 2, p)
            for i in range(n_r):
                if i != rk and Awork[i][c]:
                    f = (Awork[i][c] * inv) % p
                    Awork[i] = [(x - f * y) % p for x, y in zip(Awork[i], Awork[rk])]
                    comb_[i] = [(x - f * y) % p for x, y in zip(comb_[i], comb_[rk])]
            rk += 1
        if rk == n_r:
            break
        # rows rk..n_r-1 of Awork are zero: comb_[rk] is a dependency mod p
        coeffs = comb_[rk]
        newrow = [0] * ncol
        for j, cj in enumerate(coeffs):
            if cj:
                newrow = [x + cj * y for x, y in zip(newrow, rows[j])]
        assert all(x % p == 0 for x in newrow)
        newrow = [x // p for x in newrow]
        # replace a row j with coeff unit by newrow (keeps the Z_(p)-span saturated progressively)
        jrep = next(j for j, cj in enumerate(coeffs) if cj % p)
        rows[jrep] = newrow
    # rows now a Z_(p)-basis of the saturated lattice; but we scaled by p^Lmax: the true lattice is L, and
    # saturation of the scaled Q-space is the same lattice (saturation is scale invariant).
    return any(r[kstar_idx] % p for r in rows)

def run(d, T, m, profile, p_idx, s, p, rng, maxk=40, verbose=False):
    """profile: list of sign vectors of length d+1."""
    n = d + 1
    bs = []
    for beta in profile:
        b = [(beta[x] - beta[n - 1]) // 2 for x in range(d)]
        bs.append(b)
    pts = []
    for k in itertools.product(range(-T, T + 1), repeat=d):
        ok = True
        for qi, b in enumerate(bs):
            if qi == p_idx:
                continue
            if abs(sum(kx * bx for kx, bx in zip(k, b))) > s:
                ok = False; break
        if ok:
            pts.append(k)
    bp = bs[p_idx]
    lev = {k: sum(kx * bx for kx, bx in zip(k, bp)) for k in pts}
    mon = monos(d, m)
    E = len(mon)
    res = []
    cand = pts[:]
    rng.shuffle(cand)
    for ks in cand[:maxk]:
        lstar = lev[ks]
        S = [ks] + [k for k in pts if lev[k] > lstar]
        nlev = max(lev[k] for k in S) - lstar
        if nlev < m:
            continue
        Ev = [[binom_val(k, a, T) for a in mon] for k in S]
        r1 = rank_mod(Ev, (1 << 61) - 1)
        r0 = rank_mod(Ev[1:], (1 << 61) - 1)
        char0_sep = (r1 == r0 + 1)
        padic_sep = None
        if not char0_sep:
            levels = [lev[k] for k in S]
            padic_sep = padic_separator(S, 0, levels, lstar, mon, T, p)
        res.append((ks, lstar, nlev, len(S), char0_sep, padic_sep))
    return pts, E, res

if __name__ == '__main__':
    rng = random.Random(1)
    out = []
    configs = [
        # (d, T, m, r_narrow, s, p)
        (2, 4, 4, 2, 1, 3),
        (2, 5, 5, 2, 1, 3),
        (3, 3, 4, 2, 1, 3),
        (3, 3, 5, 3, 1, 5),
        (3, 4, 5, 3, 1, 3),
    ]
    for (d, T, m, rn, s, p) in configs:
        n = d + 1
        for trial in range(2):
            prof = [[rng.choice([-1, 1]) for _ in range(n)] for _ in range(rn + 1)]
            pts, E, res = run(d, T, m, prof, 0, s, p, rng, maxk=25)
            n0 = sum(1 for r in res if r[4])
            npad = sum(1 for r in res if (not r[4]) and r[5])
            ngood = sum(1 for r in res if (not r[4]) and (r[5] is False))
            line = (f"d={d} T={T} m={m} narrow={rn} s={s} p={p} prof={prof} |A|={len(pts)} E={E} tested={len(res)} "
                    f"char0_sep={n0} padic_sep_only={npad} CERTIFICATE(no sep)={ngood}")
            print(line); sys.stdout.flush()
            out.append(line)
            for r in res[:6]:
                print("   ", r)
