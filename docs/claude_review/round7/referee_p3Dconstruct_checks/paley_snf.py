# Referee check: Paley-I 0/1 matrix, F_2 rank and Z-surjectivity of row subsets
# (Consequence 4 of Theorem 1 needs rows R of atilde surjective mod m_q/2).
import sys, random
from fractions import Fraction

def chi(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1

def paley(q0):
    # rows: index 0 = infinity, 1..q0 = x in F_q0 ; cols: 0 = constant, 1..q0 = labels a
    M = q0 + 1
    H = [[0]*M for _ in range(M)]
    for x in range(M):
        for a in range(M):
            if x == 0 or a == 0:
                H[x][a] = 1
            else:
                xx, aa = x-1, a-1
                H[x][a] = -((1 if xx == aa else 0) + chi(aa - xx, q0))
    return H

def is_hadamard(H):
    M = len(H)
    for i in range(M):
        for j in range(M):
            s = sum(H[i][k]*H[j][k] for k in range(M))
            if s != (M if i == j else 0): return False
    return True

def rank_mod(rows, p):
    A = [r[:] for r in rows]
    n = len(A[0]) if A else 0
    rk = 0
    for c in range(n):
        piv = None
        for i in range(rk, len(A)):
            if A[i][c] % p:
                piv = i; break
        if piv is None: continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][c], -1, p)
        A[rk] = [(v*inv) % p for v in A[rk]]
        for i in range(len(A)):
            if i != rk and A[i][c] % p:
                f = A[i][c]
                A[i] = [(A[i][k] - f*A[rk][k]) % p for k in range(n)]
        rk += 1
    return rk

def det_divisor_gcd(rows):
    # surjectivity over Z of R x n matrix: compute HNF-like reduction of the column lattice:
    # surjective iff the gcd of maximal minors is 1. Use integer row-echelon on transpose.
    import math
    R = len(rows); n = len(rows[0])
    # columns as vectors in Z^R; reduce lattice basis by integer Gaussian elimination (Hermite)
    cols = [[rows[i][j] for i in range(R)] for j in range(n)]
    basis = []
    piv_row = 0
    A = cols
    # we compute the Hermite normal form of the R x n matrix via row ops on transpose
    T = [c[:] for c in A]  # n x R
    r = 0
    diag = []
    for c in range(R):
        # gcd-eliminate column c among rows r..n-1 of T
        while True:
            nz = [i for i in range(r, len(T)) if T[i][c] != 0]
            if not nz: break
            i0 = min(nz, key=lambda i: abs(T[i][c]))
            T[r], T[i0] = T[i0], T[r]
            done = True
            for i in range(r+1, len(T)):
                if T[i][c]:
                    f = T[i][c] // T[r][c]
                    T[i] = [T[i][k] - f*T[r][k] for k in range(R)]
                    if T[i][c]: done = False
            if done: break
        if r < len(T) and T[r][c] != 0:
            diag.append(abs(T[r][c])); r += 1
        else:
            diag.append(0)
    prod = 1
    for d in diag: prod *= d
    return prod  # = index of column lattice in Z^R (0 if rank < R)

if __name__ == '__main__':
    random.seed(1)
    for q0 in [43, 47, 59, 67, 71, 79, 83]:
        H = paley(q0)
        assert is_hadamard(H)
        M = q0 + 1
        P = [[(1 + H[x][a])//2 for a in range(1, M)] for x in range(M)]  # M x q0, 0/1 label columns
        r2 = rank_mod(P, 2); r3 = rank_mod(P, 3)
        # greedy: largest row subset R independent mod 2 (random order), then index over Z
        best = None
        for trial in range(20):
            order = list(range(M)); random.shuffle(order)
            R = []
            for x in order:
                if rank_mod([P[y] for y in R+[x]], 2) == len(R)+1:
                    R.append(x)
            idx = det_divisor_gcd([P[y] for y in R])
            if best is None or len(R) > best[0] or (len(R) == best[0] and idx == 1):
                best = (len(R), idx)
        print(f"q0={q0} (q0 mod 8 = {q0%8}) M={M}: rank_F2(P)={r2}, rank_F3(P)={r3}; "
              f"max F2-independent row subset size={best[0]}, Z-index of its column lattice={best[1]}")
