"""Transfer of two_thirds/proof.md Thm 5.2 to F_2[t] with the norm of F_4(t)/F_2(t).
Primes P of F_2[t] of degree e: odd e inert (residue norm fibre 2^e+1), even e split (fibre 2^e-1);
no ramified prime. log p -> deg P = e. Chord bound deg(z_i-z_j) <= D/2 + c, W = deg N = 2D.
Inequality: F2*(M) <= M*D/4 + binom(M,2) c, i.e. W >= 8 F2*(M)/M - 4(M-1)c.
"""
from math import comb, log2
import sys

def E(n, q):
    if n <= 0:
        return 0
    b, r = divmod(n, q)
    return q * comb(b, 2) + r * b

def nirr(e):
    # number of monic irreducibles of degree e over F_2 (Moebius)
    def mu(n):
        res, k, x = 1, 2, n
        while k * k <= x:
            if x % k == 0:
                x //= k
                if x % k == 0:
                    return 0
                res = -res
            k += 1
        if x > 1:
            res = -res
        return res
    return sum(mu(dd) * 2 ** (e // dd) for dd in range(1, e + 1) if e % dd == 0) // e

def depth_sum(n, q0, base):
    tot, q = 0, q0
    while q < n:
        tot += E(n, q)
        q *= base
    return tot

def psi8(M, q0, base):
    w = [8 * depth_sum(n, q0, base) for n in range(M + 1)]
    f = [None] * (M + 1)
    for S in range(1, M):
        f[S] = w[S] + (2 * S - M) ** 2
    for S in range(1, M):
        for S2 in range(S + 1, M):
            v = f[S] + w[S2 - S] + (2 * S2 - M) ** 2
            if v < f[S2]:
                f[S2] = v
    return min(f[S] + w[M - S] for S in range(1, M))

def Fstar2(M):
    tot8 = 0  # 8 * F*
    e = 1
    while True:
        q = 2 ** e
        if e % 2 == 1:
            if q + 1 >= M:
                break_inert = True
            else:
                tot8 += 8 * e * nirr(e) * depth_sum(M, q + 1, q)
        else:
            if q - 1 < M:
                free8 = 8 * depth_sum(M, q - 1, q)
                g8 = min(free8, psi8(M, q - 1, q))
                tot8 += e * nirr(e) * g8
        if q - 1 >= M:
            break
        e += 1
    return tot8 / 8

if __name__ == '__main__':
    print(" M   F2*(M)   W_min=8F*/M (c=0)   3 M log2 M   2m(m-1)(log2 m -1) [Part V deg N]  ratio")
    for M in [4, 8, 16, 32, 64, 128, 256]:
        F = Fstar2(M)
        Wmin = 8 * F / M
        r = int(round(log2(M)))
        partV = 2 * M * (M - 1) * (r - 1)
        print(f"{M:4d} {F:10.1f} {Wmin:12.1f} {3*M*log2(M):12.1f} {partV:12d} {partV/Wmin if Wmin>0 else float('inf'):10.2f}")
