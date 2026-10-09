"""Exact check of Lemma 1.2 (validity of the residue axiom (R) for actual lattice points) and of Prop 1.1's
example (the uniform Q_M strengthening of outside.md Thm A.3 fails for actual points).

(R): for points z_x on x^2+y^2 = N (N odd, no prime 3 mod 4, gcd of the points 1), a zero-sum integer c,
Z_+ = prod_{c_x>0} z_x^{c_x}, Z_- = prod_{c_x<0} z_x^{-c_x}, v = sum c_x a_x, and a unit u, let D_u(c) be the
product of the odd prime powers q^a <= Q0 (q not dividing N) with Z_+ = u Z_- mod q^a.  Then
    |Z_+ - u Z_-|^2  >=  2 D_u(c)^2 N^(n/2) / Norm(A_v),         Norm(A_v) = e^{w(v)},
which is |sin(theta - psi_u)| >= D_u e^{-w/2}/sqrt 2.  Checked in integers on random characters of
random point sets of several circles."""
import random, itertools
from gauss5 import *

def gdivides(a, z):
    """does Gaussian a divide z?"""
    n = gnorm(a)
    w = gmul(z, gconj(a))
    return w[0] % n == 0 and w[1] % n == 0

def gdiv(z, a):
    n = gnorm(a); w = gmul(z, gconj(a))
    return (w[0] // n, w[1] // n)

def points(N):
    pts = []
    from math import isqrt
    for x in range(-isqrt(N), isqrt(N) + 1):
        y2 = N - x * x
        y = isqrt(y2)
        if y * y == y2:
            pts.append((x, y))
            if y:
                pts.append((x, -y))
    return pts

random.seed(7)
Q0 = 60
pps = [(qq, a, m) for (qq, a, m) in odd_prime_powers_upto(Q0)]
units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
total = 0; tight = None
for primes in [(5, 13, 17, 29), (5, 13, 17, 29, 37), (13, 17, 29, 41, 53), (5, 13, 17, 29, 37, 41)]:
    N = 1
    for p in primes:
        N *= p
    pis = [two_squares(p) for p in primes]
    pis = [(a, b) for (a, b) in pis]
    pts = points(N)
    # primitive points only (gcd 1): not divisible by any rational prime p (i.e. not by both pi and conj pi)
    prim = []
    for z in pts:
        ok = True
        for p, pi in zip(primes, pis):
            if gdivides(pi, z) and gdivides(gconj(pi), z):
                ok = False
        if ok:
            prim.append(z)
    def expo(z):
        e = []
        for pi in pis:
            k = 0; w = z
            while gdivides(pi, w):
                w = gdiv(w, pi); k += 1
            e.append(k)
        return e
    good_pps = [(qq, a, m) for (qq, a, m) in pps if N % qq]
    for _ in range(400):
        M = random.randint(2, 6)
        Z = random.sample(prim, M)
        A = [expo(z) for z in Z]
        c = [0] * M
        while not any(c) or sum(c) != 0:
            c = [random.randint(-2, 2) for _ in range(M)]
            c[-1] -= sum(c)
        n = sum(abs(t) for t in c)
        Zp = (1, 0); Zm = (1, 0)
        for z, t in zip(Z, c):
            if t > 0:
                Zp = gmul(Zp, gpow(z, t))
            elif t < 0:
                Zm = gmul(Zm, gpow(z, -t))
        v = [sum(t * a[j] for t, a in zip(c, A)) for j in range(len(primes))]
        if not any(v):
            continue
        normA = 1
        for p, vj in zip(primes, v):
            normA *= p ** abs(vj)
        for u in units:
            D = 1
            for qq in set(t[0] for t in good_pps):
                best = 0
                for (q2, a, m) in good_pps:
                    if q2 != qq:
                        continue
                    d = (Zp[0] - (u[0] * Zm[0] - u[1] * Zm[1]), Zp[1] - (u[0] * Zm[1] + u[1] * Zm[0]))
                    if d[0] % m == 0 and d[1] % m == 0:
                        best = max(best, a)
                D *= qq ** best
            d = (Zp[0] - (u[0] * Zm[0] - u[1] * Zm[1]), Zp[1] - (u[0] * Zm[1] + u[1] * Zm[0]))
            lhs = gnorm(d) * normA
            rhs = 2 * D * D * N ** (n // 2) if n % 2 == 0 else None
            assert n % 2 == 0
            assert lhs >= rhs, (N, Z, c, u, D)
            ratio = lhs / rhs
            if D > 1 and (tight is None or ratio < tight[0]):
                tight = (ratio, N, D)
            total += 1
print("(R) validity: %d (character, unit) instances on 4 circles, all satisfy |Z+ - u Z-|^2 Norm(A_v) >= 2 D_u^2 N^(n/2);"
      " tightest ratio with D>1: %.3f (N=%d, D=%d)" % (total, tight[0], tight[1], tight[2]))

# Prop 1.1 example: the uniform Q_M strengthening (outside.md Thm A.3, property 3) fails for actual points
z1, z2 = (2, 1), (1, 2)          # on x^2 + y^2 = 5, chord sqrt 2 <= C sqrt R with C = 1.0
M = 2
from math import gcd
QM = 1
for k in range(1, 2 * M + 1):
    QM = QM * k // gcd(QM, k)
# pair character: A_v = 2+i (norm 5): property 3 requires dist(arg A_v, (pi/4)Z) >= arcsin(Q_M / sqrt 5)
print("Prop 1.1 example: points 2+i, 1+2i on x^2+y^2=5; Q_2 = lcm(1..4) = %d; property 3 needs arcsin(%d/sqrt 5) "
      "= arcsin(%.2f) -- impossible (> 1). So the uniform strengthening is not a valid input." % (QM, QM, QM / 5 ** 0.5))
