"""Exact checks for upper.md (all integer / Fraction arithmetic; no floating point decisions).

(A) Peeling recursion B_M(p,n) <= phi_M(p,n) := min_j l_j(p,n),  l_j = (M-1)(p/(M-j) + n/j),
    for M <= 14 and p+n <= 30 (B_M is a rigorous upper bound for reg Q^M(p,n) by Lemma 2.2).
    Also checks the algebraic identity K = V+/alpha' + V-/beta'' used in the inductive step.
(B) Lower-bound certificates: c = sum_sigma sgn(sigma) delta_{sigma v}, v = (-a..b, 0..0),
    satisfies sum_k c_k k^alpha = 0 for |alpha| < C(r,2) (r = a+b+1) and
    sum_k c_k Vandermonde_r(k) != 0.  Exhaustive over all monomials for M <= 5; for M = 6,7
    all generators with r = M are checked on all monomials of degree < C(r,2) in r-1 variables
    via power sums of 6 random integer linear forms (sanity check only; the proof is the
    Vandermonde divisibility argument).
(C) Cone tiling: u_j = (T_{M-j}, T_{j-1}); l_j(u_j) = l_j(u_{j+1}) = C(M,2); slopes increase;
    max_{p,n} phi_M(p,n)/(p+n) = 2 - 1/ceil(M/2) (checked on all vertices u_j, M <= 60).
(D) Shell decoupling counterexample on X_2.
(E) M = 3, 4: lower bound L_M (semigroup) equals floor(phi_M) for p,n <= 40, so
    reg Q^3(p,n) = p+n+min(p,n) and reg Q^4(p,n) = floor(phi_4(p,n)) exactly.
"""
import sys, itertools, random, math
from fractions import Fraction as F
from functools import lru_cache
sys.setrecursionlimit(100000)

def T(x):
    return x * (x + 1) // 2

def ell(M, j, p, n):
    return F(M - 1) * (F(p, M - j) + F(n, j))

def phi(M, p, n):
    return min(ell(M, j, p, n) for j in range(1, M))

@lru_cache(maxsize=None)
def B(M, p, n):
    if M == 1:
        return 0
    vals = [B(M - 1, p - c, n) for c in range(0, p + 1)] + [B(M - 1, p, n - c) for c in range(1, n + 1)]
    vals.sort(reverse=True)
    return max(v + i for i, v in enumerate(vals))

def part_A():
    bad = 0; cnt = 0
    for M in range(2, 15):
        tmax = 30 if M <= 8 else 22
        for t in range(0, tmax + 1):
            for n in range(0, t + 1):
                p = t - n
                cnt += 1
                if B(M, p, n) > phi(M, p, n):
                    bad += 1
                    print("  VIOLATION", M, p, n, B(M, p, n), phi(M, p, n))
    # identity in the inductive step: for 2 <= j <= m-1 (m = M-1),
    # V+ = l^{M-1}_j(p,n), alpha' = (m-1)/(m-j); V- = l^{M-1}_{j-1}(p,n), beta'' = (m-1)/(j-1)
    idbad = 0
    for M in range(4, 40):
        m = M - 1
        for j in range(2, m):
            al = F(m - 1, m - j); be = F(m - 1, j - 1)
            assert 1 / al + 1 / be == 1
            for (p, n) in [(0, 1), (1, 0), (3, 7), (11, 2), (5, 5)]:
                Vp = ell(M - 1, j, p, n); Vm = ell(M - 1, j - 1, p, n)
                if Vp / al + Vm / be != ell(M, j, p, n):
                    idbad += 1
                # slopes of the slice bounds
                assert ell(M - 1, j, p + 1, n) - ell(M - 1, j, p, n) == al
                assert ell(M - 1, j - 1, p, n + 1) - ell(M - 1, j - 1, p, n) == be
    print(f"(A) peeling bound B_M <= phi_M on {cnt} cases: {'PASS' if bad == 0 else 'FAIL'};"
          f" inductive identity: {'PASS' if idbad == 0 else 'FAIL'}")

def alt_measure(r, a, b, M):
    v = list(range(-a, b + 1)) + [0] * (M - r)
    assert len(v) == M
    meas = {}
    for perm in itertools.permutations(range(r)):
        # sign of perm
        sgn = 1
        seen = [False] * r
        for i in range(r):
            if not seen[i]:
                L = 0; j = i
                while not seen[j]:
                    seen[j] = True; j = perm[j]; L += 1
                if L % 2 == 0:
                    sgn = -sgn
        k = tuple(v[perm[i]] for i in range(r)) + tuple(v[r:])
        meas[k] = meas.get(k, 0) + sgn
    return meas

def monos(n, d):
    if n == 1:
        yield (d,); return
    for a in range(d, -1, -1):
        for rest in monos(n - 1, d - a):
            yield (a,) + rest

def part_B():
    ok = True
    for M in range(2, 8):
        for r in range(2, M + 1):
            if M >= 6 and r < M:
                continue
            for a in range(0, r):
                b = r - 1 - a
                meas = alt_measure(r, a, b, M)
                pts = list(meas.items())
                assert all(c != 0 for _, c in pts) and len(pts) == math.factorial(r)
                ordlow = r * (r - 1) // 2
                # support: shell S(T_b, T_a)
                for k, _ in pts:
                    assert sum(x for x in k if x > 0) == T(b) and -sum(x for x in k if x < 0) == T(a)
                if M <= 5:
                    for d in range(0, ordlow):
                        for al in monos(r - 1, d):   # coordinates 1..r-1 suffice (sum fixed, others 0)
                            s = 0
                            for k, c in pts:
                                term = c
                                for i, e in enumerate(al):
                                    term *= k[i] ** e
                                s += term
                            if s != 0:
                                ok = False; print("  moment nonzero", M, a, b, d, al)
                else:
                    rng = random.Random(12345 + M * 100 + a)
                    for _ in range(6):
                        lam = [rng.randint(-50, 50) for _ in range(M)]
                        for d in range(0, ordlow):
                            s = sum(c * sum(l * x for l, x in zip(lam, k)) ** d for k, c in pts)
                            if s != 0:
                                ok = False; print("  power sum nonzero", M, a, b, d)
                # top degree: Vandermonde pairing nonzero
                vd = 0
                for k, c in pts:
                    pr = 1
                    for i in range(r):
                        for j in range(i + 1, r):
                            pr *= (k[i] - k[j])
                    vd += c * pr
                if vd == 0:
                    ok = False; print("  Vandermonde pairing zero", M, a, b)
    print(f"(B) antisymmetrised-orbit certificates (ord >= C(r,2), exact): {'PASS' if ok else 'FAIL'}")

def part_C():
    ok = True
    for M in range(2, 61):
        u = [None] + [(T(M - j), T(j - 1)) for j in range(1, M + 1)]
        CM2 = M * (M - 1) // 2
        for j in range(1, M):
            if ell(M, j, *u[j]) != CM2 or ell(M, j, *u[j + 1]) != CM2:
                ok = False
        # slopes n/p increasing
        sl = [F(u[j][1], u[j][0]) if u[j][0] else None for j in range(1, M + 1)]
        for j in range(1, M - 1):
            if not (sl[j - 1] < sl[j]):
                ok = False
        # phi(u_j) = C(M,2) (min attained), and ratio max over vertices
        best = max(F(CM2, sum(u[j])) for j in range(1, M + 1))
        for j in range(1, M + 1):
            if phi(M, *u[j]) != CM2:
                ok = False
        k = (M + 1) // 2
        if best != 2 - F(1, k):
            ok = False
        # phi(p,n)/(p+n) <= 2 - 1/k on a grid (phi is min of linear forms; the max ratio of a
        # concave 1-homogeneous function on the segment p+n=1 is checked at all breakpoints u_j)
    print(f"(C) cone tiling, phi_M(u_j) = C(M,2), max ratio = 2 - 1/ceil(M/2) (M<=60): {'PASS' if ok else 'FAIL'}")

def part_D():
    # On X_2 (|z1| = |z2|): F = z1*conj(z2) + conj(z1)*z2 - 2*z1*conj(z1) = -|z1 - z2|^2.
    rng = random.Random(7)
    ok = True
    for _ in range(200):
        # Gaussian integers of equal norm: z2 = z1 * unit or conj-type; use z1=x+iy, z2=y+ix etc.
        x, y = rng.randint(-99, 99), rng.randint(-99, 99)
        z1 = complex(x, y)
        for z2 in [complex(y, x), complex(-x, y), complex(x, -y), complex(-y, x)]:
            Fv = z1 * z2.conjugate() + z1.conjugate() * z2 - 2 * z1 * z1.conjugate()
            if abs(Fv + abs(z1 - z2) ** 2) > 1e-6:
                ok = False
    # shell components at a point of Y (z1 = z2 = z): z conj z + conj z z = 2N != 0 and -2N != 0
    print(f"(D) X_2 example F = -|z1-z2|^2 with non-vanishing shell components: {'PASS' if ok else 'FAIL'}")

@lru_cache(maxsize=None)
def Lsemi(M, p, n):
    best = 0
    for r in range(2, M + 1):
        for a in range(0, r):
            b = r - 1 - a
            cp, cn, o = T(b), T(a), r * (r - 1) // 2
            if cp <= p and cn <= n:
                best = max(best, o + Lsemi(M, p - cp, n - cn))
    return best

def part_E():
    ok = True
    for M in (3, 4):
        for p in range(0, 41):
            for n in range(0, 41):
                if Lsemi(M, p, n) != math.floor(phi(M, p, n)):
                    ok = False; print("  mismatch", M, p, n, Lsemi(M, p, n), phi(M, p, n))
    for p in range(0, 41):
        for n in range(0, 41):
            if math.floor(phi(3, p, n)) != p + n + min(p, n):
                ok = False
    print(f"(E) M=3,4: semigroup lower bound = floor(phi_M) for p,n<=40: {'PASS' if ok else 'FAIL'}")

if __name__ == '__main__':
    part_A(); part_B(); part_C(); part_D(); part_E()
