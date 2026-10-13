# Referee's own construction of slab-box models (independent of the author's code).
# Exact arithmetic: python-flint fmpz_mat / nmod_mat.
import itertools, sys
import flint

def sylvester_cuts(k):
    """Sylvester-Hadamard of order n=2^k: points x in F_2^k (index = integer), cuts S_a={x: a.x=0}, a!=0."""
    n = 2**k
    cuts = []
    for a in range(1, n):
        S = frozenset(x for x in range(n) if bin(a & x).count('1') % 2 == 0)
        cuts.append(S)
    return n, cuts

def slab_box(n, cuts, T, L, base=None):
    """A = {k' in Z^{n-1}: |k_x|<=T (x != base), |k(S)|<=L for all cuts}, k_base = -sum.
    Returns list of full vectors k in Z^n (sum 0)."""
    if base is None:
        base = n - 1
    others = [x for x in range(n) if x != base]
    A = []
    for kp in itertools.product(range(-T, T + 1), repeat=n - 1):
        k = [0] * n
        for x, v in zip(others, kp):
            k[x] = v
        k[base] = -sum(kp)
        ok = True
        for S in cuts:
            s = sum(k[x] for x in S)
            if abs(s) > L:
                ok = False
                break
        if ok:
            A.append(tuple(k))
    return A, others

def reduced_monomials(d, deg_lt, maxexp):
    """exponent vectors alpha in N^d with alpha_x <= maxexp and |alpha| < deg_lt"""
    out = []
    for al in itertools.product(range(maxexp + 1), repeat=d):
        if sum(al) < deg_lt:
            out.append(al)
    return out

def moment_matrix(A, others, m, T):
    # On A, coordinate k_x (x in others) takes values in [-T,T]; prod_{|s|<=T}(k_x - s) vanishes,
    # so monomials with alpha_x <= 2T span Poly_{<m}|_A (reduction only lowers degree).
    d = len(others)
    mons = reduced_monomials(d, m, 2 * T)
    rows = []
    for al in mons:
        row = []
        for k in A:
            v = 1
            for x, e in zip(others, al):
                if e:
                    v *= k[x] ** e
            row.append(v)
        rows.append(row)
    return rows, mons

def rank_mod(rows, p=2**61 - 1):
    if not rows:
        return 0
    M = flint.nmod_mat(len(rows), len(rows[0]), [v % p for r in rows for v in r], p)
    return M.rank()

def rank_exact(rows):
    if not rows:
        return 0
    M = flint.fmpz_mat(rows)
    return M.rank()

def rho_and_dims(A, others, T, mmax, exact=False):
    """dim Poly_{<m}|_A for m=1..mmax; rho = max m with dim V(A,m)>0."""
    res = {}
    for m in range(1, mmax + 1):
        rows, _ = moment_matrix(A, others, m, T)
        r = rank_exact(rows) if exact else rank_mod(rows)
        res[m] = r
        if r == len(A):
            break
    return res

if __name__ == '__main__':
    pass
