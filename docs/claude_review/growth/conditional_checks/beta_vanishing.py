"""Average vanishing order along the small diagonal Y of sections of O(N) on X_M.

X_M = {u_1 v_1 = ... = u_M v_M} (over C; the circle form is a twist).  Degree-N
sections of O(N)|X_M have the basis t^k * prod u_j^{e_j} (e_j>0) prod v_j^{-e_j}
(e_j<0), 2k+|e|_1 = N.  Near Y (u_j = u w_j, v_j = v / w_j, w_1 = 1) such a monomial
restricts to u^a v^b * prod_{j>=2} w_j^{e_j}; grouping by s = sum_j e_j, the order of
vanishing along Y of a section is min over s of the order at w=1 of a Laurent
polynomial supported on
   A_s = { e' in Z^{M-1} : |s - sum e'| + |e'|_1 <= N },   s = N mod 2.
Order >= m  <=>  coefficients annihilate all polynomials of degree < m on A_s.
Hence  h^0(NH - mE) = sum_s (|A_s| - h_{A_s}(m-1)), and the total order
  sum_{m>=1} h^0(NH-mE) = sum_s T(A_s),  T(A) = sum_{d>=0} (|A| - h_A(d)).
beta_N := sum_s T(A_s) / (N * sum_s |A_s|)  ->  beta(H,E) = int vol(H-tE)dt / vol(H).
h_A(d) = rank of degree<=d monomials evaluated on A (computed mod two primes).
"""
import itertools, sys
import numpy as np

def points(n, N, s):
    pts = []
    rng = range(-N, N + 1)
    for e in itertools.product(rng, repeat=n):
        if abs(s - sum(e)) + sum(abs(x) for x in e) <= N:
            pts.append(e)
    return pts

def monomials_of_degree(n, d):
    # exponent vectors of total degree d
    if n == 1:
        yield (d,)
        return
    for a in range(d + 1):
        for rest in monomials_of_degree(n - 1, d - a):
            yield (a,) + rest

def hilbert_T(pts, p):
    A = np.array(pts, dtype=np.int64) % p  # coordinates mod p
    size = len(pts)
    n = A.shape[1]
    # basis in reduced form: rows of B, pivot columns
    B = np.zeros((0, size), dtype=np.int64)
    piv = []
    T = 0
    d = 0
    hist = []
    # powers cache
    powcache = {}
    def colpow(j, a):
        key = (j, a)
        if key not in powcache:
            if a == 0:
                powcache[key] = np.ones(size, dtype=np.int64)
            else:
                powcache[key] = (colpow(j, a - 1) * A[:, j]) % p
        return powcache[key]
    rank = 0
    while rank < size:
        for mono in monomials_of_degree(n, d):
            v = np.ones(size, dtype=np.int64)
            for j, a in enumerate(mono):
                if a:
                    v = (v * colpow(j, a)) % p
            if rank:
                coeff = v[piv]  # since B is in reduced form
                v = (v - (coeff @ B) % p) % p
            nz = np.nonzero(v)[0]
            if len(nz) == 0:
                continue
            c = nz[0]
            inv = pow(int(v[c]), p - 2, p)
            v = (v * inv) % p
            # eliminate column c from existing rows
            if rank:
                col = B[:, c].copy()
                B = (B - np.outer(col, v) % p) % p
            B = np.vstack([B, v])
            piv.append(c)
            rank += 1
            if rank == size:
                break
        hist.append(rank)
        T += size - rank
        d += 1
    return T, hist

def beta(n, N, p=1000003):
    totT = 0
    totA = 0
    for s in range(-N, N + 1):
        if (s - N) % 2:
            continue
        if s < 0:
            continue
        pts = points(n, N, s)
        T, hist = hilbert_T(pts, p)
        w = 1 if s == 0 else 2
        totT += w * T
        totA += w * len(pts)
    return totT, totA, totT / (N * totA)

if __name__ == '__main__':
    n = int(sys.argv[1]); Ns = [int(x) for x in sys.argv[2:]]
    for N in Ns:
        T1, A1, b1 = beta(n, N, 1000003)
        T2, A2, b2 = beta(n, N, 999983)
        print(f"M={n+1} N={N} h0={A1} sumT={T1} (check {T2}) beta_N={b1:.5f}", flush=True)
