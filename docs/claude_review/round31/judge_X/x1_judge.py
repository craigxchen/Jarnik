#!/usr/bin/env python3
"""Judge's own check of X.1 (order lemma, mechanism) and of the X.5 row P=100 (least cutoff 439).
Exact integer arithmetic only. Imports nothing from Astra or the referees."""
import random, math, sys
from fractions import Fraction as F

def isprime(n):
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True

def gauss_factor(p):
    # a + bi with a^2+b^2 = p
    a = 1
    while True:
        b2 = p - a * a
        b = math.isqrt(b2)
        if b * b == b2: return (a, b)
        a += 1

def gmul(u, v): return (u[0]*v[0]-u[1]*v[1], u[0]*v[1]+u[1]*v[0])
def gconj(u): return (u[0], -u[1])
def gpow(u, n):
    r = (1, 0)
    for _ in range(n): r = gmul(r, u)
    return r

# ---------- (1) order lemma brute force ----------
def order_lemma():
    bad = 0; tests = 0
    for N in range(2, 300):
        for k in range(0, 6):
            n = 2 ** k
            for u in range(N):  # additive cyclic group Z/N
                d = N // math.gcd(u, N)           # ord(u)
                dn = d // math.gcd(d, n)          # ord(n u)
                if dn == 4:
                    tests += 1
                    if d != 4 * n: bad += 1
    print("order lemma: tests", tests, "failures", bad)

# ---------- (2) mechanism on random lifts ----------
def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def mechanism(S, trials, H, rng):
    T = max(v2(l - 1) for l in S)
    n = 2 ** (T - 1)
    fails = 0; cnt = 0
    for _ in range(trials):
        p = rng.choice(S)
        pi = gauss_factor(p)
        if rng.random() < .5: pi = gconj(pi)
        # random gamma with Norm <= H
        while True:
            g = (rng.randint(-int(H**.5), int(H**.5)), rng.randint(-int(H**.5), int(H**.5)))
            if 0 < g[0]**2 + g[1]**2 <= H: break
        Ng = g[0]**2 + g[1]**2
        num = gmul(pi, gmul(g, g))   # r = pi gamma^2 / Norm gamma
        gg = math.gcd(math.gcd(num[0], num[1]), Ng)
        alpha = (num[0] // gg, num[1] // gg); b = Ng // gg
        if b < 0: b = -b; alpha = (-alpha[0], -alpha[1])
        J = alpha[0]**2 + alpha[1]**2
        assert J == p * b * b
        if b == 1: continue
        cnt += 1
        x, y = gpow(alpha, n)
        ok = True
        for D in (x - y, x + y):
            if D == 0 or D % 2 == 0: ok = False
            if math.gcd(D, J) != 1: ok = False
            for l in S:
                if D % l == 0: ok = False
            if abs(D) > math.isqrt(2 * J ** n) + 1: ok = False
        # sum of squares identity: (x-y)^2+(x+y)^2 = 2 J^n, so RHS sum = J^n/p^n = b^{2n}
        if (x - y) ** 2 + (x + y) ** 2 != 2 * J ** n: ok = False
        if not ok: fails += 1
    print(f"mechanism: |S|={len(S)} T={T} n={n} H={H} nonintegral lifts {cnt} failures {fails}")

# ---------- (3) table row P=100 ----------
def oddext_L(X, S):
    L = 1
    for q in range(3, X + 1, 2):
        if isprime(q) and q not in S:
            a = q
            while a * q <= X: a *= q
            L *= a
    return L

def killed_by_n(alpha, p, n, L):
    x, y = gpow(alpha, n)
    D = [y, x - y, x, x + y]
    m = [math.gcd(d, L) for d in D]
    pn = p ** n
    # constraints: cos2t in [-1+2 m2^2/pn, 1-2 m0^2/pn], sin2t in [-1+m3^2/pn, 1-m1^2/pn]
    c_lo = F(-pn + 2 * m[2] ** 2, pn); c_hi = F(pn - 2 * m[0] ** 2, pn)
    s_lo = F(-pn + m[3] ** 2, pn); s_hi = F(pn - m[1] ** 2, pn)
    if c_lo > c_hi or s_lo > s_hi: return True
    # does rectangle meet unit circle? min dist^2 <= 1 <= max dist^2
    def clamp0(lo, hi):
        if lo <= 0 <= hi: return F(0)
        return min(abs(lo), abs(hi))
    mn = clamp0(c_lo, c_hi) ** 2 + clamp0(s_lo, s_hi) ** 2
    mx = max(c_lo ** 2, c_hi ** 2) + max(s_lo ** 2, s_hi ** 2)
    return not (mn <= 1 <= mx)

def survives_phi(alpha, p, L, phi_w=(2, 1)):
    """phi = arg(2+i): check all n until p^n >= 5^n L^2 (then automatic)."""
    Nw = phi_w[0]**2 + phi_w[1]**2
    n = 0
    while True:
        n += 1
        x, y = gpow(alpha, n)
        D = [y, x - y, x, x + y]
        m = [math.gcd(d, L) for d in D]
        X_, Y_ = gpow(phi_w, n)
        E = [Y_, Y_ - X_, X_, X_ + Y_]   # sin(n phi - s pi/4) * nu_s^{-1}... up to sign
        # inequality: |E_s| nu_s / Nw^{n/2} >= nu_s m_s / p^{n/2}  <=>  E_s^2 p^n >= m_s^2 Nw^n
        for s in range(4):
            if E[s] ** 2 * p ** n < m[s] ** 2 * Nw ** n: return False, n
        if p ** n >= Nw ** n * L * L: return True, n

def table_row():
    P = 100
    S = [p for p in range(P, 2 * P + 1) if isprime(p) and p % 4 == 1]
    cases = []
    for p in S:
        pp = gauss_factor(p)
        for l in S:
            pl = gauss_factor(l)
            for plc in (pl, gconj(pl)):
                alpha = gmul(pp, gmul(plc, plc))
                g = math.gcd(alpha[0], alpha[1])
                if g != 1: continue        # non-primitive case l = p
                cases.append((p, l, alpha))
    print("P=100 labels", len(S), "cases", len(cases))
    res = {}
    for X in (437, 438, 439):
        L = oddext_L(X, set(S))
        alive = []
        for (p, l, alpha) in cases:
            if not any(killed_by_n(alpha, p, n, L) for n in range(1, 41)):
                alive.append((p, l, alpha))
        res[X] = alive
        print(f"X={X}: cases not killed by any single n<=40: {len(alive)}", alive[:5])
    # explicit survivor at 438
    L = oddext_L(438, set(S))
    for (p, l, alpha) in res[438]:
        for conjflag in (False, True):
            a = gconj(alpha) if conjflag else alpha
            ok, nmax = survives_phi(a, p, L)
            print("  survivor check", p, l, a, "phi=arg(2+i) passes all n:", ok, "checked to n =", nmax)

if __name__ == "__main__":
    order_lemma()
    rng = random.Random(1)
    S100 = [p for p in range(100, 201) if isprime(p) and p % 4 == 1]
    S58 = [p for p in range(1000, 2001) if isprime(p) and p % 8 == 5]
    mechanism(S100, 3000, 400, rng)
    mechanism(S58, 3000, 4000, rng)
    table_row()
