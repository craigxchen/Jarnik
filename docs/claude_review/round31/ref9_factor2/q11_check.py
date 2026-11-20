"""Referee's own exact check of the IX.3 q=11 computation (factor2 block).
Exact Gaussian-rational arithmetic in Q(i); no import of Astra's scripts.
Usage: python3 q11_check.py [mine|flint]   (choice of integral right inverse U)
"""
import sys
from itertools import combinations
from math import isqrt, gcd
from sympy import primerange, factorint, isprime

# ---------- Gaussian integers as (a,b) ----------
def gm(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gc(x): return (x[0], -x[1])
def gn(x): return x[0]*x[0]+x[1]*x[1]
def gpow(x, e):
    r = (1, 0)
    for _ in range(e): r = gm(r, x)
    return r
IU = (0, 1)
def ipow(e): return [(1,0),(0,1),(-1,0),(0,-1)][e % 4]

# Gaussian rationals as (num, den) with num, den Gaussian integers
def qmul(a, b): return (gm(a[0], b[0]), gm(a[1], b[1]))
def qpow(a, e):
    if e >= 0: return (gpow(a[0], e), gpow(a[1], e))
    return (gpow(a[1], -e), gpow(a[0], -e))
def qconj(a): return (gc(a[0]), gc(a[1]))
def qeq(a, b): return gm(a[0], b[1]) == gm(b[0], a[1])
def qinv(a): return (a[1], a[0])
def red(a, mod):
    """reduce Gaussian rational modulo integer mod (den must be a unit mod mod)"""
    n, d = a
    nd = gn(d)
    inv = pow(nd % mod, -1, mod)
    x = gm(n, gc(d))
    return ((x[0]*inv) % mod, (x[1]*inv) % mod)
def rmul(x, y, mod): z = gm(x, y); return (z[0] % mod, z[1] % mod)

# ---------- integer matrices ----------
def matmul(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

def right_inverse(D):
    """Unimodular column reduction: D Q = [L|0], L lower triangular unimodular; U = Q[:, :r] L^{-1}."""
    r, c = len(D), len(D[0])
    M = [row[:] for row in D]
    Q = [[int(i == j) for j in range(c)] for i in range(c)]
    def colop(j, k, t):  # col_j -= t col_k
        for row in M: row[j] -= t*row[k]
        for row in Q: row[j] -= t*row[k]
    def swap(j, k):
        for row in M: row[j], row[k] = row[k], row[j]
        for row in Q: row[j], row[k] = row[k], row[j]
    for p in range(r):
        while True:
            nz = [j for j in range(p, c) if M[p][j] != 0]
            if not nz: raise ValueError("rank deficient")
            jmin = min(nz, key=lambda j: abs(M[p][j]))
            swap(p, jmin)
            done = True
            for j in range(p+1, c):
                if M[p][j]:
                    colop(j, p, M[p][j] // M[p][p]);
                    if M[p][j]: done = False
            if done: break
        if abs(M[p][p]) != 1: raise ValueError("Smith invariant not one")
    L = [row[:r] for row in M]
    # invert lower triangular unimodular L
    Linv = [[0]*r for _ in range(r)]
    for col in range(r):
        e = [int(i == col) for i in range(r)]
        x = [0]*r
        for i in range(r):
            s = e[i] - sum(L[i][k]*x[k] for k in range(i))
            assert s % L[i][i] == 0
            x[i] = s // L[i][i]
        for i in range(r): Linv[i][col] = x[i]
    U = matmul([row[:r] for row in Q], Linv)
    return U

def main(which):
    q = 11
    labels = list(range(1, q))
    sq = {(i*i) % q for i in range(1, q)}
    chi = lambda t: 0 if t % q == 0 else (1 if t % q in sq else -1)
    K = lambda i, j: chi(j - i)
    marks = [i for i in labels if K(i, 0) == 1]
    m, nl = len(marks), len(labels)
    A = [[1 + K(i, j) for j in labels] for i in marks]
    assert all(set(col) >= {0, 2} for col in zip(*A)), "endpoint coverage"
    for x, i in enumerate(marks):   # odd exponent exactly at own label
        assert [j for j in labels if A[x][j-1] % 2] == [i]
    D = [[A[x][j] - A[0][j] for j in range(nl)] for x in range(1, m)]
    H = [[int(labels[j] == marks[x]) - int(labels[j] == marks[0]) for j in range(nl)] for x in range(1, m)]
    assert all(sum(abs(t) for t in row) == q + 1 for row in D)
    assert all((a - b) % 2 == 0 for r1, r2 in zip(D, H) for a, b in zip(r1, r2))
    if which == 'flint':
        from flint import fmpz_mat
        Dm = fmpz_mat(D)
        herm, T = Dm.transpose().hnf(transform=True)
        U = [[int(T[i, j]) for i in range(m-1)] for j in range(nl)]
    else:
        U = right_inverse(D)
    DU = matmul(D, U)
    assert DU == [[int(i == j) for j in range(m-1)] for i in range(m-1)]
    HD2 = [[(H[i][j] - D[i][j]) // 2 for j in range(nl)] for i in range(m-1)]
    V = matmul(U, HD2)
    F = [[int(i == j) + 2*V[i][j] for j in range(nl)] for i in range(nl)]
    assert matmul(D, F) == H
    print("q=11 marks", marks, "U shape", len(U), "x", len(U[0]), "max|U| =", max(abs(t) for r in U for t in r),
          "max|V| =", max(abs(t) for r in V for t in r), "max|F| =", max(abs(t) for r in F for t in r))
    # primes as in Astra's script (recovered from it): first 10 primes = 1 mod 4 from 101; factor with least x>=1
    primes = [p for p in primerange(101, 500) if p % 4 == 1][:nl]
    pis = []
    for p in primes:
        x = next(x for x in range(1, isqrt(p)+1) if isqrt(p-x*x)**2 == p-x*x)
        pis.append((x, isqrt(p - x*x)))
    gam = [(pi, gc(pi)) for pi in pis]   # gamma_l = pi_l / conj(pi_l)
    winding = [2*x*x + 4 for x in range(m)]
    dk = [(winding[x] - winding[0]) // 2 for x in range(1, m)]
    u = [sum(U[j][x]*dk[x] for x in range(m-1)) for j in range(nl)]
    assert [sum(D[x][j]*u[j] for j in range(nl)) for x in range(m-1)] == dk

    def witnesses(skip=None):
        rr = []
        for j in range(nl):
            r = ((ipow(u[j]) if True else None), (1, 0))
            r = qmul(r, (pis[j], (1, 0)))
            for l in range(nl):
                if l == skip: continue
                r = qmul(r, qpow(gam[l], V[j][l]))
            rr.append(r)
        return rr
    rr = witnesses()
    for j, r in enumerate(rr):     # exact norm and exact ratio identity (IX.6)
        assert gn(r[0]) == primes[j]*gn(r[1])
        rat = qmul(r, qinv(qconj(r)))
        target = ((-1)**(u[j] % 2), 1)
        tg = ((target[0], 0), (1, 0))
        for l in range(nl): tg = qmul(tg, qpow(gam[l], F[j][l]))
        assert qeq(rat, tg)
    # exact point marks in Q(i) and exact pair ratio identity
    def marksQ(rr, skip=None):
        out = []
        for x in range(m):
            z = (ipow(-winding[x]), (1, 0))
            for j in range(nl):
                if j == skip: continue
                z = qmul(z, qpow(rr[j], A[x][j]))
                z = qmul(z, qpow(qconj(rr[j]), 2 - A[x][j]))
            out.append(z)
        return out
    pts = marksQ(rr)
    for x, y in combinations(range(m), 2):
        lhs = qmul(pts[x], qinv(pts[y]))
        rhs = qmul(gam[marks[x]-1], qinv(gam[marks[y]-1]))
        assert qeq(lhs, rhs)
    print("exact Q(i) identity rho_x/rho_y = gamma_x/gamma_y holds for all", m*(m-1)//2, "pairs")
    # modular checks as in the claim
    ext = 0
    for ell in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        for a in (1, 2, 3):
            b = ell**a
            P = [red(z, b) for z in pts]
            G = [red(g, b) for g in gam]
            for x, y in combinations(range(m), 2):
                assert rmul(P[x], G[marks[y]-1], b) == rmul(P[y], G[marks[x]-1], b)
                ext += 1
    same = 0
    for l, p in enumerate(primes):
        rl = witnesses(skip=l)
        for j in range(nl):
            if j != l: assert gn(rl[j][0]) == primes[j]*gn(rl[j][1])
        lv = marksQ(rl, skip=l)
        for a in (1, 2, 3):
            b = p**a
            P = [red(z, b) for z in lv]
            G = [red(g, b) if gi != l else None for gi, g in enumerate(gam)]
            for x, y in combinations(range(m), 2):
                if A[x][l] == A[y][l]:
                    assert l + 1 not in (marks[x], marks[y])
                    assert rmul(P[x], G[marks[y]-1], b) == rmul(P[y], G[marks[x]-1], b)
                    # exact as well
                    assert qeq(qmul(lv[x], qinv(lv[y])), qmul(gam[marks[x]-1], qinv(gam[marks[y]-1])))
                    same += 1
    print("exterior pair identities:", ext, " same-level identities:", same)
    # J_j and alpha_j for each row of F
    Js = []
    for j in range(nl):
        J = 1
        for l in range(nl): J *= primes[l]**abs(F[j][l])
        Js.append(J)
    print("rows with J_j > p_j:", sum(J > p for J, p in zip(Js, primes)), "of", nl,
          "; log10(J_j/p_j) =", [round(len(str(J//p)), 0) for J, p in zip(Js, primes)])
    # singleton test: ratio r_j/conj r_j = i^e mod ell, ell^2 > 2 p_j (s odd) or ell^2 > p_j (s even)
    found = []
    for ell in primerange(3, 2000):
        if ell in primes: continue
        for j, r in enumerate(rr):
            rat = red(qmul(r, qinv(qconj(r))), ell)
            for e in range(4):
                ie = ipow(e); ie = (ie[0] % ell, ie[1] % ell)
                if rat == ie:
                    nu2 = 2 if e % 2 else 1     # 1/nu_s^2
                    if ell*ell > nu2*primes[j]:
                        found.append((ell, primes[j], e))
    print("first singleton e-all failures (ell, p_j, s):", found[:6])
    print("contains (23,173,3):", (23, 173, 3) in found)
    # monomial n e_j, s = 0, n = 0 mod 4 : exterior part of Im(alpha^n) vs p_j^{n/2}
    parent = set(primes)
    for j in range(nl):
        al = (1, 0)
        for l in range(nl):
            f = F[j][l]
            al = gm(al, gpow(pis[l], f) if f > 0 else gpow(gc(pis[l]), -f))
        assert gn(al) == Js[j]
        # check ratio
        first = None
        for n in range(4, 41, 4):
            Y = gpow(al, n)[1]
            assert Y != 0
            Y = abs(Y)
            while Y % 2 == 0: Y //= 2
            for p in parent:
                while Y % p == 0: Y //= p
            if Y > primes[j]**(n//2):
                first = n; break
        if j < 3 or first is None:
            print("  row", j, "p_j", primes[j], "first n=0 mod 4 (<=40) with exterior part of Im(alpha^n) > p_j^(n/2):", first)
    return

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'mine')
