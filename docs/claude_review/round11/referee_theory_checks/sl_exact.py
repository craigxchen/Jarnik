"""Exact SL-full tester (referee's own code; no author code imported).

Given a finite S in Z^d (points k'), an order m, an integer level function l on S and a prime p, decide for
every k* in S whether some c in V(S,m) (rational vectors annihilating all polynomials of degree < m on S)
has a unique minimiser of v_p(c_k) + l(k) at k*  (property UM_p at k*).
Method (theory.md Prop 3.4(i), re-derived): Lambda = {(c_k p^{l(k)-lmin})} cap Z_p^S is saturated; UM_p at k*
iff e_{k*} lies in R = Lambda mod p.  Exact rational nullspace, then p-saturation by the standard loop
(replace a row b_j by (sum lambda_i b_i)/p while the rows are dependent mod p), all in Python integers.
"""
from fractions import Fraction
from math import gcd
from itertools import combinations_with_replacement

def monomial_exps(d, m):
    out = []
    for deg in range(m):
        for combo in combinations_with_replacement(range(d), deg):
            a = [0] * d
            for i in combo:
                a[i] += 1
            out.append(tuple(a))
    return out

def moment_matrix(S, m):
    d = len(S[0])
    exps = monomial_exps(d, m)
    M = []
    for a in exps:
        row = []
        for k in S:
            v = 1
            for x in range(d):
                if a[x]:
                    v *= k[x] ** a[x]
            row.append(v)
        M.append(row)
    return M

def nullspace_Q(M, ncols):
    """rational basis of {c : M c = 0}, returned as integer primitive vectors."""
    A = [[Fraction(v) for v in r] for r in M]
    piv_cols = []
    rk = 0
    for c in range(ncols):
        piv = None
        for i in range(rk, len(A)):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = 1 / A[rk][c]
        A[rk] = [x * inv for x in A[rk]]
        for i in range(len(A)):
            if i != rk and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[rk])]
        piv_cols.append(c)
        rk += 1
    free = [c for c in range(ncols) if c not in set(piv_cols)]
    basis = []
    for fcol in free:
        v = [Fraction(0)] * ncols
        v[fcol] = Fraction(1)
        for i, pc in enumerate(piv_cols):
            v[pc] = -A[i][fcol]
        den = 1
        for x in v:
            den = den * x.denominator // gcd(den, x.denominator)
        iv = [int(x * den) for x in v]
        g = 0
        for x in iv:
            g = gcd(g, abs(x))
        basis.append([x // g for x in iv])
    return basis

def rref_mod_p(rows, p):
    A = [[x % p for x in r] for r in rows]
    ncols = len(A[0]) if A else 0
    rk = 0
    pivs = []
    for c in range(ncols):
        piv = None
        for i in range(rk, len(A)):
            if A[i][c]:
                piv = i
                break
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][c], p - 2, p)
        A[rk] = [(x * inv) % p for x in A[rk]]
        for i in range(len(A)):
            if i != rk and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[rk])]
        pivs.append(c)
        rk += 1
    return A[:rk], pivs

def dependency_mod_p(rows, p):
    """return lambda (list) with sum lambda_i rows_i == 0 mod p, lambda != 0 mod p, or None."""
    r = len(rows)
    ncols = len(rows[0])
    # augment with identity to track combinations
    A = [[x % p for x in rows[i]] + [1 if j == i else 0 for j in range(r)] for i in range(r)]
    rk = 0
    for c in range(ncols):
        piv = None
        for i in range(rk, r):
            if A[i][c]:
                piv = i
                break
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][c], p - 2, p)
        A[rk] = [(x * inv) % p for x in A[rk]]
        for i in range(r):
            if i != rk and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[rk])]
        rk += 1
    if rk == r:
        return None
    return A[rk][ncols:]   # a zero row of the left block: its tracking part is a dependency

def p_saturate(B, p):
    B = [list(r) for r in B]
    while True:
        lam = dependency_mod_p(B, p)
        if lam is None:
            return B
        j = next(i for i, x in enumerate(lam) if x % p)
        comb = [sum(lam[i] * B[i][c] for i in range(len(B))) for c in range(len(B[0]))]
        assert all(x % p == 0 for x in comb)
        B[j] = [x // p for x in comb]

def um_points(S, m, lvals, p):
    """return (dimV, list of indices k* with UM_p achievable in V(S,m))."""
    M = moment_matrix(S, m)
    V = nullspace_Q(M, len(S))
    if not V:
        return 0, []
    lmin = min(lvals)
    D = [p ** (l - lmin) for l in lvals]
    B = [[v[k] * D[k] for k in range(len(S))] for v in V]
    B = p_saturate(B, p)
    R, pivs = rref_mod_p(B, p)
    good = []
    pivset = set(pivs)
    for row, pc in zip(R, pivs):
        if all(row[c] == 0 for c in range(len(S)) if c != pc):
            good.append(pc)
    return len(V), sorted(good)

def vp(x, p):
    from fractions import Fraction as F
    x = F(x)
    if x == 0:
        return None
    v = 0
    a, b = x.numerator, x.denominator
    while a % p == 0:
        a //= p; v += 1
    while b % p == 0:
        b //= p; v -= 1
    return v

def certificate(S, m, lvals, p, kstar):
    """Construct an explicit c in V(S,m) with UM_p at kstar (if one exists) and verify it from scratch:
    M c = 0 exactly, and v_p(c_k)+l(k) has a unique minimiser at kstar.  Returns (ok, c)."""
    from fractions import Fraction as F
    M = moment_matrix(S, m)
    V = nullspace_Q(M, len(S))
    lmin = min(lvals)
    D = [p ** (l - lmin) for l in lvals]
    B = p_saturate([[v[k] * D[k] for k in range(len(S))] for v in V], p)
    # solve lambda B == e_kstar (mod p): transpose system  B^T lambda = e
    nb = len(B); ns = len(S)
    A = [[B[i][k] % p for i in range(nb)] + [1 if k == kstar else 0] for k in range(ns)]
    rk = 0; pivc = []
    for c in range(nb):
        piv = next((i for i in range(rk, ns) if A[i][c]), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][c], p - 2, p)
        A[rk] = [(x * inv) % p for x in A[rk]]
        for i in range(ns):
            if i != rk and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[rk])]
        pivc.append(c); rk += 1
    if any(A[i][nb] for i in range(rk, ns)):
        return False, None
    lam = [0] * nb
    for i, c in enumerate(pivc):
        lam[c] = A[i][nb]
    cp = [sum(lam[i] * B[i][k] for i in range(nb)) for k in range(ns)]
    c = [F(cp[k], D[k]) for k in range(ns)]
    # verify from scratch
    for row in M:
        assert sum(F(a) * b for a, b in zip(row, c)) == 0
    vals = [(vp(c[k], p) + lvals[k], k) for k in range(ns) if c[k] != 0]
    mn = min(v for v, _ in vals)
    arg = [k for v, k in vals if v == mn]
    return (arg == [kstar]), c
