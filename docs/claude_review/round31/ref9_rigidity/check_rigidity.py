"""Referee checks for IX.7 (rigidity). Exact integer arithmetic only (except
final float comparisons of cosines, which are double-checked by exact
rational inequalities where they matter).
"""
from math import gcd, lcm, isqrt
from fractions import Fraction
import random, math
from sympy import factorint, isprime, primerange

def mul(a, b, m=None):
    r = (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
    return (r[0] % m, r[1] % m) if m else r
def gpow(a, n, m=None):
    r = (1, 0)
    for _ in range(n):
        r = mul(r, a, m)
    return r
def conj(a): return (a[0], -a[1])
def norm(a): return a[0]**2+a[1]**2
def inv_mod(a, m):
    nn = norm(a) % m
    t = pow(nn, -1, m)
    c = conj(a)
    return (c[0]*t % m, c[1]*t % m)

# ---------------------------------------------------------------- (3) control
z = (2, 1)
z6 = gpow(z, 6)
assert z6 == (-117, 44), z6
assert factorint(117) == {3: 2, 13: 1}
for Q in (9, 13):
    ratio = mul(z, inv_mod(conj(z), Q), Q)
    assert gpow(ratio, 6, Q) == ((-1) % Q, 0)
    # maximality at 27 (for 9) and 169 (for 13)
for Q in (27, 169):
    ratio = mul(z, inv_mod(conj(z), Q), Q)
    assert gpow(ratio, 6, Q) != ((-1) % Q, 0)
# cos(2 beta) = Re((2+i)^2)/5 = 3/5 exactly; need 3/5 >= 117/125 for e-all
assert Fraction(3, 5) < Fraction(117, 125)
print("control: (2+i)^6 =", z6, "; exterior unit -1 collision part 9*13 = 117; 3/5 < 117/125: e-all FAILS at X>=13")

# exact identity |Re((2+i)^6)|/5^3 = |cos(6 beta)| ; and cos(6*beta/3)=cos(2beta)=3/5
beta = math.atan2(1, 2)
assert abs(abs(math.cos(6*beta/3)) - 0.6) < 1e-12

# --------------------------------------------- (b) clean real parts, general alpha
def order_mod(a, ell):
    w = (a[0] % ell, a[1] % ell); o = 1
    while w != (1 % ell, 0):
        w = mul(w, a, ell); o += 1
        if o > 10**6: raise RuntimeError
    return o

def clean_check(alpha, E, K=25):
    J = norm(alpha)
    L = 4
    for ell in E:
        if J % ell:
            L = lcm(L, order_mod(alpha, ell))
    for k in range(1, K+1):
        re = gpow(alpha, L*k)[0]
        for ell in E:
            assert re % ell != 0, (alpha, ell, k)
        if J % 2 == 1:
            assert re % 2 == 1
    return L

# alpha = (2+i)*(3+2i)^2 : one side of each pair, J = 5*13^2
a1 = mul((2, 1), gpow((3, 2), 2))
for E in [(2, 3, 5, 7, 13), (2, 5, 13, 17, 29, 3, 11)]:
    L = clean_check(a1, E)
    print("clean real part OK for alpha =", a1, "E =", E, "L =", L)

# ------------------------------- (c) fake rational lift fails e-all at unit -1
# r = (2+i)(3+2i)/(3-2i): norm 5, r/conj r = alpha/conj alpha with alpha=(2+i)(3+2i)^2
# Data: residues r mod q^a; angle phi = arg r (the "natural" angle for this lift).
# e-all at v = n e_1, s=2:  |cos(n phi)| >= m_{n,2} * 5^{-n/2}, where m = prod over
# exterior q^a (q odd, q not 5, q not 13 [exceptional: support of r]) with
# (r/conj r)^n = -1 mod q^a.  Since (r/conj r)^n = (alpha/conj alpha)^n,
# collision at q^a iff q^a | Re(alpha^n) (alpha unit mod q).
def exterior_part(x, E):
    x = abs(x)
    for ell in E:
        while x % ell == 0:
            x //= ell
    return x

alpha = a1
E = (2, 5, 13)
phi = math.atan2(alpha[1], alpha[0])  # arg alpha = arg r^2... see note
# arg(r) : r = (2+i)(3+2i)^2/13 ; arg r = arg alpha exactly
r_num = alpha  # r = alpha/13 since (3+2i)/(3-2i) = (3+2i)^2/13
first_fail = None
for n in range(1, 80):
    an = gpow(alpha, n)
    m = exterior_part(an[0], E)
    # |cos(n phi)| = |Re alpha^n| / J^{n/2}
    J = norm(alpha)
    # test m * 5^{-n/2} <= |Re alpha^n| / J^{n/2}  <=>  m^2 * J^n <= Re^2 * 5^n
    ok = m*m*J**n <= an[0]**2 * 5**n
    if not ok and first_fail is None:
        first_fail = (n, m, an[0])
print("fake lift r=(2+i)(3+2i)/(3-2i), phi=arg r: first e-all(unit -1) failure at n, m, Re:", first_fail)
assert first_fail is not None

# The same lift with the genuine angle phi = arg(2+i) also fails (theorem says any angle fails).
phi2 = math.atan2(1, 2)
ff = None
for n in range(1, 200):
    an = gpow(alpha, n)
    m = exterior_part(an[0], E)
    lhs = abs(math.cos(n*phi2)); rhs = m * 5**(-n/2)
    if rhs > lhs*(1+1e-9):
        ff = (n, m); break
print("same lift with phi = arg(2+i): first failure at n, m =", ff)

# ------------------------------------------------- genuine passes (sanity)
alpha = (2, 1)
for n in range(1, 120):
    an = gpow(alpha, n)
    m = exterior_part(an[0], (2, 5))
    assert m*m <= an[0]**2   # |cos n beta| = |Re|/5^{n/2} >= m/5^{n/2}
print("genuine 2+i with phi=beta passes rays n<=119 (trivially m <= |Re|)")

# ------------------- weak approximation at a finite cutoff: any (d)-data lifts
# For p=5, pi=2+i, and odd prime powers Q=q^{A_q} <= X (q != 5), choose random
# units r(q) of Z[i]/Q with Norm r(q) = 5 mod Q, then find gamma in Z[i] with
# pi*gamma/conj(gamma) = r(q) mod Q for all q, by CRT.
def random_data(X, seed):
    rnd = random.Random(seed)
    data = {}
    for q in primerange(3, X+1):
        if q == 5: continue
        A = 1
        while q**(A+1) <= X: A += 1
        Q = q**A
        while True:
            a = (rnd.randrange(Q), rnd.randrange(Q))
            if gcd(norm(a), q) != 1: continue
            # scale to norm 5: need 5/Norm(a) a square mod Q
            t = 5 * pow(norm(a), -1, Q) % Q
            # find s with s^2 = t mod Q (brute force, small Q)
            s = next((s for s in range(Q) if s*s % Q == t), None)
            if s is None: continue
            data[q] = (Q, ((a[0]*s) % Q, (a[1]*s) % Q))
            break
    return data

def lift(data, pi=(2, 1)):
    # want gamma with gamma/conj(gamma) = r(q)/pi mod Q.  Norm(r/pi) = 1, and
    # by Prop 2.1(2) t = g/conj g with g a unit; brute force g mod Q, then CRT.
    rems = []
    for q, (Q, rq) in data.items():
        t = mul(rq, inv_mod(pi, Q), Q)
        g = None
        for x in range(Q):
            for y in range(Q):
                if gcd(x*x+y*y, q) != 1: continue
                if mul((x, y), inv_mod((x, -y), Q), Q) == t:
                    g = (x, y); break
            if g: break
        assert g is not None
        rems.append((Q, g))
    # CRT componentwise
    from sympy.ntheory.modular import crt
    mods = [Q for Q, _ in rems]
    gx = crt(mods, [g[0] for _, g in rems])[0]
    gy = crt(mods, [g[1] for _, g in rems])[0]
    gam = (int(gx), int(gy))
    for q, (Q, rq) in data.items():
        val = mul(mul(pi, gam, Q), inv_mod(conj(gam), Q), Q)
        assert val == rq, (q, val, rq)
    return gam

for seed in range(5):
    d = random_data(40, seed)
    g = lift(d)
    print("cutoff 40, seed", seed, ": data lifted by r = (2+i)*gamma/conj(gamma), norm(gamma) has",
          len(str(norm(g))), "digits")
print("ALL CHECKS DONE")
