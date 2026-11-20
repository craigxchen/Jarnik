#!/usr/bin/env python3
"""Judge check of the IX.3 witnesses and the same-level deletion against sharp.md Def. 1.1 / lower.md Def. 1.1(e).

Own code: Paley matrix, own unimodular column reduction for an integral right inverse U of D,
V = U(H-D)/2, F = I + 2V, witnesses r_j = i^{u_j} pi_j prod_l gamma_l^{V_jl} in Z[i]/Q,
level witnesses at p_l with gamma_l omitted.  Checks, for several q and random S, windings with even
differences:
 (1) DU = I, DF = H, V integral;
 (2) at exterior Q = l^a (l odd prime, not a label, a<=3): Norm r_j = p_j, every pair ratio
     rho_x/rho_y = gamma_x/gamma_y, random zero-sum characters prod rho^c = prod gamma^c;
 (3) at every parent prime power p_l^a (a<=3): Norm r_j^(l) = p_j (j != l); level residues have norm N/p_l^2;
     for every SAME-LEVEL pair, rho_x(p_l)/rho_y(p_l) = gamma_x/gamma_y, and the record's collision
     exponent max{a<=3 : rho_x = rho_y mod p_l^a} equals min(3, v_{p_l}(Im(pi_x conj pi_y)));
 (4) the owner l is the unique mark at level 1 at label l;
 (5) negative control: odd winding difference breaks the exterior identity.
"""
import random, sys
from sympy import primerange, isprime
from math import isqrt

def leg(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1

def col_reduce_right_inverse(D):
    """D: list of rows (m x n), full row rank, Smith invariants 1. Return integral U (n x m) with DU=I."""
    m, n = len(D), len(D[0])
    A = [row[:] for row in D]
    Q = [[1 if i == j else 0 for j in range(n)] for i in range(n)]  # column ops tracked: A = D Q
    def colop_sub(dst, src, k):  # col dst -= k*col src
        for r in range(m): A[r][dst] -= k*A[r][src]
        for r in range(n): Q[r][dst] -= k*Q[r][src]
    def swap(c1, c2):
        for r in range(m): A[r][c1], A[r][c2] = A[r][c2], A[r][c1]
        for r in range(n): Q[r][c1], Q[r][c2] = Q[r][c2], Q[r][c1]
    def neg(c):
        for r in range(m): A[r][c] = -A[r][c]
        for r in range(n): Q[r][c] = -Q[r][c]
    for r in range(m):
        # gcd-reduce row r over columns r..n-1 into column r
        while True:
            nz = [c for c in range(r, n) if A[r][c] != 0]
            if len(nz) == 1 and nz[0] == r: break
            if not nz: raise ValueError("rank")
            piv = min(nz, key=lambda c: abs(A[r][c]))
            if piv != r: swap(piv, r)
            for c in range(r+1, n):
                if A[r][c]: colop_sub(c, r, A[r][c]//A[r][r])
        if A[r][r] < 0: neg(r)
        assert A[r][r] == 1, ("pivot", A[r][r])
    # A = [L | 0], L lower triangular unit diagonal; U = Q[:, :m] L^{-1}
    L = [row[:m] for row in A]
    # invert L (unit lower triangular) over Z
    Linv = [[0]*m for _ in range(m)]
    for j in range(m):
        e = [1 if i == j else 0 for i in range(m)]
        x = [0]*m
        for i in range(m):
            x[i] = e[i] - sum(L[i][k]*x[k] for k in range(i))
        for i in range(m): Linv[i][j] = x[i]
    U = [[sum(Q[r][k]*Linv[k][c] for k in range(m)) for c in range(m)] for r in range(n)]
    return U

def matmul(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

# Gaussian arithmetic mod Q
def mul(a, b, Q): return ((a[0]*b[0]-a[1]*b[1]) % Q, (a[0]*b[1]+a[1]*b[0]) % Q)
def conj(a, Q): return (a[0] % Q, (-a[1]) % Q)
def inv(a, Q):
    n = (a[0]*a[0]+a[1]*a[1]) % Q
    ni = pow(n, -1, Q)
    return ((a[0]*ni) % Q, (-a[1]*ni) % Q)
def pw(a, k, Q):
    if k < 0: a, k = inv(a, Q), -k
    r = (1, 0)
    while k:
        if k & 1: r = mul(r, a, Q)
        a = mul(a, a, Q); k >>= 1
    return r
def ipow(k, Q): return [(1, 0), (0, 1), (Q-1, 0), (0, Q-1)][k % 4]

def split_pi(p):
    for x in range(1, isqrt(p)+1):
        y2 = p-x*x; y = isqrt(y2)
        if y*y == y2: return (x, y)

def run(q, rho, seed, odd_control=False):
    rnd = random.Random(seed)
    I = list(range(1, q))
    A = [i for i in I if leg(0 - i, q) == 1]       # K_{i0} = chi(0-i) = 1
    S = [i for i in A if rnd.random() < rho]
    if len(S) < 3: return None
    m = len(S)
    K = lambda i, j: leg(j - i, q)
    a = {i: [1 + K(i, j) for j in I] for i in S}
    # coverage
    for c in range(q-1):
        vals = {a[i][c] for i in S}
        if not (0 in vals and 2 in vals): return ('nocover',)
    i0 = S[0]
    D = [[a[i][c] - a[i0][c] for c in range(q-1)] for i in S[1:]]
    H = [[(1 if I[c] == i else 0) - (1 if I[c] == i0 else 0) for c in range(q-1)] for i in S[1:]]
    U = col_reduce_right_inverse(D)
    DU = matmul(D, U)
    assert all(DU[r][c] == (r == c) for r in range(m-1) for c in range(m-1))
    HmD = [[H[r][c]-D[r][c] for c in range(q-1)] for r in range(m-1)]
    assert all(x % 2 == 0 for row in HmD for x in row)
    V = matmul(U, [[x//2 for x in row] for row in HmD])     # (q-1) x (q-1)
    F = [[(1 if r == c else 0) + 2*V[r][c] for c in range(q-1)] for r in range(q-1)]
    DF = matmul(D, F)
    assert DF == H
    # primes: distinct split primes randomly assigned
    pr = [p for p in primerange(101, 5000) if p % 4 == 1]
    rnd.shuffle(pr); P = pr[:q-1]
    pis = [split_pi(p) for p in P]
    # windings with even differences (nonconstant)
    base = rnd.randrange(4)
    k = {i: base + 2*rnd.randrange(-3, 4) for i in S}
    if odd_control: k[S[1]] += 1
    kd = [(k[i]-k[i0])//2 for i in S[1:]]
    u = [sum(U[r][c]*kd[c] for c in range(m-1)) for r in range(q-1)]
    res = dict(q=q, m=m, ext_pair=0, ext_char=0, lev=0, lev_exp=0, fail=0, maxU=max(abs(x) for row in U for x in row))
    labels = set(P)
    def gamma(l, Q): return mul(pis[l], inv(conj(pis[l], Q), Q), Q)
    # exterior moduli
    ext = [l for l in primerange(3, 60) if l not in labels]
    for l in ext:
        for e in (1, 2, 3):
            Q = l**e
            g = [gamma(t, Q) for t in range(q-1)]
            r = []
            for j in range(q-1):
                x = mul(ipow(u[j], Q), pis[j], Q)
                for t in range(q-1):
                    if V[j][t]: x = mul(x, pw(g[t], V[j][t], Q), Q)
                assert mul(x, conj(x, Q), Q) == (P[j] % Q, 0)
                r.append(x)
            rhox = {}
            for i in S:
                x = ipow(-k[i], Q)
                for j in range(q-1):
                    x = mul(x, mul(pw(r[j], a[i][j], Q), pw(conj(r[j], Q), 2-a[i][j], Q), Q), Q)
                rhox[i] = x
            for x_ in S:
                for y_ in S:
                    if x_ == y_: continue
                    lhs = mul(rhox[x_], inv(rhox[y_], Q), Q)
                    rhs = mul(g[x_-1], inv(g[y_-1], Q), Q)
                    res['ext_pair'] += 1
                    if lhs != rhs: res['fail'] += 1
            for _ in range(20):
                c = [rnd.randrange(-2, 3) for _ in S]; c[0] -= sum(c)
                lhs = (1, 0); rhs = (1, 0)
                for t, i in enumerate(S):
                    lhs = mul(lhs, pw(rhox[i], c[t], Q), Q); rhs = mul(rhs, pw(g[i-1], c[t], Q), Q)
                res['ext_char'] += 1
                if lhs != rhs: res['fail'] += 1
    if odd_control: return res
    # parent primes
    for l in range(q-1):
        pl = P[l]
        lab = I[l]
        lev1 = [i for i in S if a[i][l] == 1]
        assert lev1 == ([lab] if lab in S else []), (lev1, lab)
        for e in (1, 2, 3):
            Q = pl**e
            g = {t: gamma(t, Q) for t in range(q-1) if t != l}
            r = {}
            for j in range(q-1):
                if j == l: continue
                x = mul(ipow(u[j], Q), pis[j], Q)
                for t in range(q-1):
                    if t != l and V[j][t]: x = mul(x, pw(g[t], V[j][t], Q), Q)
                assert mul(x, conj(x, Q), Q) == (P[j] % Q, 0)
                r[j] = x
            Nl = 1
            for j in range(q-1):
                if j != l: Nl = Nl*P[j]*P[j] % Q
            rhox = {}
            for i in S:
                x = ipow(-k[i], Q)
                for j in r:
                    x = mul(x, mul(pw(r[j], a[i][j], Q), pw(conj(r[j], Q), 2-a[i][j], Q), Q), Q)
                assert mul(x, conj(x, Q), Q) == (Nl, 0)
                rhox[i] = x
            for x_ in S:
                for y_ in S:
                    if x_ == y_ or a[x_][l] != a[y_][l]: continue
                    lhs = mul(rhox[x_], inv(rhox[y_], Q), Q)
                    rhs = mul(g[x_-1], inv(g[y_-1], Q), Q)
                    res['lev'] += 1
                    if lhs != rhs: res['fail'] += 1
    # record collision exponent at p_l for same-level pairs, through cubes
    for l in range(q-1):
        pl = P[l]
        for x_ in S:
            for y_ in S:
                if x_ >= y_ or a[x_][l] != a[y_][l]: continue
                px, py = pis[x_-1], pis[y_-1]
                im = px[1]*py[0] - px[0]*py[1]   # Im(px * conj py)
                vv = 0
                while im % pl == 0: im //= pl; vv += 1
                # direct collision exponent from residues
                aa = 0
                for e in (1, 2, 3):
                    Q = pl**e
                    # recompute level residues (cheap enough for small q)
                    g = {t: gamma(t, Q) for t in range(q-1) if t != l}
                    def lr(i):
                        x = ipow(-k[i], Q)
                        for j in range(q-1):
                            if j == l: continue
                            rj = mul(ipow(u[j], Q), pis[j], Q)
                            for t in range(q-1):
                                if t != l and V[j][t]: rj = mul(rj, pw(g[t], V[j][t], Q), Q)
                            x = mul(x, mul(pw(rj, a[i][j], Q), pw(conj(rj, Q), 2-a[i][j], Q), Q), Q)
                        return x
                    if lr(x_) == lr(y_): aa = e
                    else: break
                res['lev_exp'] += 1
                if aa != min(3, vv): res['fail'] += 1
    return res

if __name__ == '__main__':
    out = []
    for q in (11, 19, 23):
        for seed in range(3):
            for rho in (1.0, 0.8):
                r = run(q, rho, seed*7+q)
                print(r, flush=True)
    print('odd control', run(19, 1.0, 5, odd_control=True))
