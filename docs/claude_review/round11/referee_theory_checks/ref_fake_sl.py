"""Referee check (exact integer arithmetic): the barrier profiles violate Hypothesis SL(n) of theory.md.

Profiles built independently from the definitions in round5/construct.md Sec 2 (Paley-I, identity skew
pairing, b unflipped copies + one flipped copy per row) and round5/exact.md Sec 4 (skew double core).
For each profile we check:
  (a) admissibility ingredients: exponent rank n-1, pair separation (equal weights) >= 1/2 - 1/(2n);
  (b) the twin structure: for a column p, a 'twin' is a column q != p with a_p - a_q supported on a set X
      of non-base points with |X| <= 2.  Then on A_p (slab of q kept narrow, |l_q| <= w/2) the form
      l_p = l_q + sum_{x in X} +-k_x has width <= w + 2T|X| <= w + 4T < m (regime: m >= 4w), so l_p takes
      <= m values on A_p and the level-sum lemma excludes UM_p: SL(P,p) fails.
      We count, for every choice of base point, the columns WITHOUT a twin ('escaping' columns).
SL(n) needs >= 8 SL columns; escaping <= 7 for some base point refutes SL(n) at that n (in the regime).
"""
from fractions import Fraction
import sys

def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1

def paley_I(q):
    # rows/cols indexed by inf (=0) and F_q (=1..q). H(inf,.) = 1, H(i,inf) = -1, H(i,j) = delta + chi(j-i)
    M = q + 1
    H = [[0] * M for _ in range(M)]
    for j in range(M):
        H[0][j] = 1
    for i in range(1, M):
        H[i][0] = -1
        for j in range(1, M):
            H[i][j] = (1 if i == j else 0) + legendre((j - 1) - (i - 1), q)
    # normalise: multiply rows i in F_q by -1  -> column 0 all ones
    for i in range(1, M):
        H[i] = [-v for v in H[i]]
    return H

def check_hadamard(H):
    M = len(H)
    for i in range(M):
        for j in range(M):
            s = sum(H[i][k] * H[j][k] for k in range(M))
            assert s == (M if i == j else 0)

def round5_profile(q, b=8):
    H = paley_I(q)
    check_hadamard(H)
    M = q + 1
    assert all(H[x][0] == 1 for x in range(M))
    # rows x; labels 1..M-1; pairing pi(i) = i for i in F_q (row i <-> label i), x_inf = row 0, a* = label 1
    a = {x: x for x in range(1, M)}
    a[0] = 1
    cols = []  # each column is a list of signs over rows
    for lab in range(1, M):
        for _ in range(b):
            cols.append([H[x][lab] for x in range(M)])
    for x in range(M):
        c = [H[y][a[x]] for y in range(M)]
        c[x] = -c[x]
        cols.append(c)
    return cols, M

def skew_double_core(q):
    # skew-Hadamard K of order M = q+1 from Paley: K = I + S, S antisymmetric; normalised so K(x,0) = 1
    M = q + 1
    S = [[0] * M for _ in range(M)]
    for j in range(1, M):
        S[0][j] = 1
        S[j][0] = -1
    for i in range(1, M):
        for j in range(1, M):
            if i != j:
                S[i][j] = legendre((j - 1) - (i - 1), q)
    K = [[(1 if i == j else 0) + S[i][j] for j in range(M)] for i in range(M)]
    # check KK^T = M I and K + K^T = 2I
    for i in range(M):
        for j in range(M):
            assert sum(K[i][k] * K[j][k] for k in range(M)) == (M if i == j else 0)
            assert K[i][j] + K[j][i] == (2 if i == j else 0)
    # normalise K(x,0)=1 via DKD
    d = [K[x][0] for x in range(M)]
    K = [[d[i] * K[i][j] * d[j] for j in range(M)] for i in range(M)]
    assert all(K[x][0] == 1 for x in range(M))
    cols = []
    for a_ in range(1, M):
        Sa = [1 if K[x][a_] == -1 else 0 for x in range(M)]
        Spa = Sa[:]
        Spa[0] ^= 1
        Spa[a_] ^= 1
        cols.append([2 * v - 1 for v in Sa])
        cols.append([2 * v - 1 for v in Spa])
    return cols, M

def rank_Q(rows):
    A = [[Fraction(v) for v in r] for r in rows]
    rk = 0
    ncol = len(A[0]) if A else 0
    for c in range(ncol):
        piv = None
        for i in range(rk, len(A)):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        for i in range(len(A)):
            if i != rk and A[i][c] != 0:
                f = A[i][c] / A[rk][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[rk])]
        rk += 1
    return rk

def analyse(name, cols, M):
    n = M
    r = len(cols)
    expo = [[(cols[p][x] + 1) // 2 for p in range(r)] for x in range(n)]  # a_x in {0,1}^r
    diffs = [[expo[x][p] - expo[n - 1][p] for p in range(r)] for x in range(n - 1)]
    rk = rank_Q(diffs)
    minsep = min(sum(1 for p in range(r) if cols[p][x] != cols[p][y]) for x in range(n) for y in range(x + 1, n))
    frac = Fraction(minsep, r)
    print(f"{name}: n={n} r={r} exponent rank={rk} (n-1={n-1}); min #columns separating a pair={minsep}, "
          f"fraction={float(frac):.4f} vs 1/2-1/(2n)={0.5 - 0.5 / n:.4f}  admissible-sep: {frac >= Fraction(1, 2) - Fraction(1, 2 * n)}")
    best = None
    for base in range(n):
        esc = 0
        for p in range(r):
            found = False
            for qq in range(r):
                if qq == p:
                    continue
                X = [x for x in range(n) if cols[p][x] != cols[qq][x]]
                if len(X) <= 2 and base not in X:
                    found = True
                    break
            if not found:
                esc += 1
        if best is None or esc < best[1]:
            best = (base, esc)
    print(f"   min over base points of #columns without a twin: {best[1]} (base point row {best[0]}); SL(n) needs >= 8")
    return best[1]

if __name__ == '__main__':
    for q in (7, 11, 19, 23):
        cols, M = round5_profile(q)
        analyse(f"round5 Paley fake q={q}", cols, M)
    for q in (7, 11, 19, 23):
        cols, M = skew_double_core(q)
        analyse(f"skew double core q={q}", cols, M)
