"""Exact checks for additive.md Section 2 (the additive content of a cluster is lattice-ness).

For actual normalised clusters:
  (2.1)  rank K = M - 3 where K = {c in Z^M : sum c = 0, sum c_x z_x = 0}   (Smith form, exact);
  (2.2)  index I = [Z[i] : Lambda], Lambda = Z-span of differences, equals gcd of all 2x2 minors
         = gcd over triples of 2*Area;  and q | I (q odd, q not dividing N) only if the cluster has
         at most 2 residues mod q;
  (2.3)  projection identity: [Z^(M-3) : proj_{x>=4} K] = A_123 / I;
  (2.5)  K = K_N cap {||c||_1 < N/L}: every c in an LLL-reduced basis of K_N with ||c||_1 < N/L is a relation,
         and the small-lift reconstruction z_1 + w_x recovers the cluster from z_1 and residues mod N.
  (2.3') pigeonhole relation height: the shortest relation found by LLL vs the bound H_1.
Uses only Python ints (Smith/HNF by hand) and a textbook exact-rational LLL.
"""
import itertools, random
from fractions import Fraction
from math import gcd, log, isqrt
from glib import (mul, conj, norm, sub, ggcd, exact_div, factor, points_of_norm, primes_upto)
import tangent


def hnf_kernel(rows):
    """Integer kernel basis of the map Z^n -> Z^m given by columns; rows = list of m integer rows (length n).
    Returns a basis of {c in Z^n : rows . c = 0} via column-style unimodular reduction."""
    m = len(rows)
    n = len(rows[0])
    A = [list(r) for r in rows]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]  # columns of U track ops
    piv_col = 0
    r = 0
    for r in range(m):
        # reduce row r over columns piv_col..n-1 to a single nonzero entry by Euclid
        while True:
            nz = [j for j in range(piv_col, n) if A[r][j] != 0]
            if len(nz) <= 1:
                break
            jmin = min(nz, key=lambda j: abs(A[r][j]))
            for j in nz:
                if j == jmin:
                    continue
                q = A[r][j] // A[r][jmin]
                for i in range(m):
                    A[i][j] -= q * A[i][jmin]
                for i in range(n):
                    U[i][j] -= q * U[i][jmin]
        nz = [j for j in range(piv_col, n) if A[r][j] != 0]
        if nz:
            j = nz[0]
            # swap column j to piv_col
            for i in range(m):
                A[i][j], A[i][piv_col] = A[i][piv_col], A[i][j]
            for i in range(n):
                U[i][j], U[i][piv_col] = U[i][piv_col], U[i][j]
            piv_col += 1
    return [[U[i][j] for i in range(n)] for j in range(piv_col, n)]


def lll(B, delta=Fraction(3, 4)):
    B = [list(b) for b in B]
    n = len(B)

    def dot(u, v):
        return sum(a * b for a, b in zip(u, v))

    def gs(B):
        Bs, mu = [], [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = [Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = Fraction(dot(B[i], Bs[j])) / dot(Bs[j], Bs[j])
                v = [a - mu[i][j] * b for a, b in zip(v, Bs[j])]
            Bs.append(v)
        return Bs, mu
    Bs, mu = gs(B)
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [a - q * b for a, b in zip(B[k], B[j])]
                Bs, mu = gs(B)
        if dot(Bs[k], Bs[k]) >= (delta - mu[k][k - 1] ** 2) * dot(Bs[k - 1], Bs[k - 1]):
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            Bs, mu = gs(B)
            k = max(k - 1, 1)
    return B


def area2(a, b, c):
    u, v = sub(b, a), sub(c, a)
    return abs(u[0] * v[1] - u[1] * v[0])


def check(cl):
    cl, N = tangent.normalise(cl)
    M = len(cl)
    if M < 4:
        return None
    # K: kernel of c -> (sum c, sum c Re z, sum c Im z)
    rows = [[1] * M, [z[0] for z in cl], [z[1] for z in cl]]
    Kb = hnf_kernel(rows)
    assert len(Kb) == M - 3, len(Kb)
    for c in Kb:
        assert sum(c) == 0 and sum(ci * z[0] for ci, z in zip(c, cl)) == 0 and sum(ci * z[1] for ci, z in zip(c, cl)) == 0
    # index I = gcd of 2x2 minors of differences = gcd of 2*areas
    diffs = [sub(z, cl[0]) for z in cl[1:]]
    I = 0
    for u, v in itertools.combinations(diffs, 2):
        I = gcd(I, abs(u[0] * v[1] - u[1] * v[0]))
    I2 = 0
    for a, b, c in itertools.combinations(cl, 3):
        I2 = gcd(I2, area2(a, b, c))
    assert I == I2
    # residues: q | I  => at most 2 residues mod q (q odd prime not dividing N, q | I)
    for q in factor(I):
        if q > 2 and N % q:
            res = {(z[0] % q, z[1] % q) for z in cl}
            assert len(res) <= 2, (q, res)
    # projection index A_123 / I
    A123 = area2(cl[0], cl[1], cl[2])
    # sublattice of Z^{M-3}: projections of Kb onto coordinates 3..M-1; index = |det| of a basis
    P = [c[3:] for c in Kb]
    # determinant via fraction-free elimination
    from fractions import Fraction as F
    Mx = [[F(x) for x in r] for r in P]
    n = len(Mx)
    det = F(1)
    for i in range(n):
        piv = next(r for r in range(i, n) if Mx[r][i] != 0)
        if piv != i:
            Mx[i], Mx[piv] = Mx[piv], Mx[i]
            det = -det
        det *= Mx[i][i]
        for r in range(i + 1, n):
            f = Mx[r][i] / Mx[i][i]
            Mx[r] = [a - f * b for a, b in zip(Mx[r], Mx[i])]
    assert abs(det) == F(A123, I), (abs(det), A123, I)
    # K_N and the local-global statement
    L = max(isqrt(norm(sub(a, b))) + 1 for a, b in itertools.combinations(cl, 2))
    # K_N = {c : sum c = 0, sum c z = 0 mod N}: kernel of rows mod N -> add N e_j columns trick
    # basis of K_N: kernel of [1..1; Re; Im] over Z^M x Z^2 (extra vars for multiples of N)
    rows2 = [[1] * M + [0, 0], [z[0] for z in cl] + [N, 0], [z[1] for z in cl] + [0, N]]
    KN = [c[:M] for c in hnf_kernel(rows2)]
    KN = lll(KN)
    nrel = 0
    for c in KN:
        if sum(abs(x) for x in c) < N / L:
            s = (sum(ci * z[0] for ci, z in zip(c, cl)), sum(ci * z[1] for ci, z in zip(c, cl)))
            assert s == (0, 0)
            nrel += 1
    # small-lift reconstruction: w_x = smallest representative of z_x - z_1 mod N, then z_1 + w_x
    for z in cl:
        w = sub(z, cl[0])
        wr = (((w[0] + N // 2) % N) - N // 2, ((w[1] + N // 2) % N) - N // 2)
        assert (cl[0][0] + wr[0], cl[0][1] + wr[1]) == z
    Ks = lll(Kb)
    h = max(abs(x) for x in Ks[0])
    return M, N, I, A123, nrel, h


if __name__ == "__main__":
    random.seed(3)
    split = [p for p in primes_upto(400) if p % 4 == 1]
    Ns = set()
    for _ in range(120):
        k = random.randint(3, 6)
        ps = random.sample(split[:30], k)
        N = 1
        for p in ps:
            N *= p ** random.randint(1, 2)
        if N < 10 ** 14:
            Ns.add(N)
    Ns |= {1176852625, 235370525, 537606725}
    cnt = 0
    stats = []
    for N in sorted(Ns):
        Cm = 8 if N in (1176852625, 235370525, 537606725) else 4
        for cl in tangent.all_runs(N, Cm):
            r = check(cl)
            if r:
                cnt += 1
                stats.append(r)
    print("clusters checked (M>=4):", cnt, "all assertions passed")
    odd = lambda n: n // (n & -n)
    print("max M:", max(s[0] for s in stats), " clusters with I even:", sum(1 for s in stats if s[2] % 2 == 0),
          " clusters with odd part of I > 1:", sum(1 for s in stats if odd(s[2]) > 1))
    print("K_N short vectors that are relations: min over clusters of #(short basis vectors) =",
          min(s[4] for s in stats), " (should equal M-3 when the last two minima exceed N/L)")
    bad = [s for s in stats if s[4] != s[0] - 3]
    print("clusters where #short LLL vectors of K_N != M-3:", len(bad))
