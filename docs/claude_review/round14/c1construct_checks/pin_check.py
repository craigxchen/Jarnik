# Checks for Prop P (norm-1 pinning in regime (1.4)).
# (a) |sin(theta - s pi/4)| = nu_s |D_s(pi)| / sqrt(p) (exact via squares: compare 2*sin^2 etc. rationally)
# (b) solution set of the one-sided conditions |sin(x - s pi/4)| >= |sin(theta - s pi/4)|, s in S,
#     x mod pi, is {theta} (b<a, S={2,3}) or {theta, pi/2-theta} (b>a, S={1,2,3}); grid evidence.
# (c) how often the p-part condition (i) fails for random prime sets.
import sys, math, random
sys.path.insert(0, '.')
from lib14 import *
from fractions import Fraction
import numpy as np

# (a) exact: sin^2(theta - s pi/4) as a rational: sin theta = b/sqrt p, cos theta = a/sqrt p.
# s=0: b^2/p ; s=1: (b-a)^2/(2p); s=2: a^2/p; s=3: (a+b)^2/(2p). nu_s^2 D_s^2/p: s even D^2/p, s odd D^2/(2p).
ok_a = 0
for p in split_primes(5, 100000):
    a, b = gauss_prime(p)
    sin2 = [Fraction(b*b, p), Fraction((b-a)**2, 2*p), Fraction(a*a, p), Fraction((a+b)**2, 2*p)]
    for s in range(4):
        D = Ds((a, b), s)
        nus2 = Fraction(1) if s % 2 == 0 else Fraction(1, 2)
        assert sin2[s] == nus2 * Fraction(D*D, p)
    # direct float sanity
    th = math.atan2(b, a)
    for s in range(4):
        assert abs(abs(math.sin(th - s*math.pi/4))**2 - float(sin2[s])) < 1e-12
    ok_a += 1
print(f"(a) identity sin^2(theta - s pi/4) = nu_s^2 D_s^2/p exact for {ok_a} split primes <= 1e5: OK")

# (b) grid evidence for the solution set
def sols(a, b, S, n=200000):
    p = a*a + b*b
    th = math.atan2(b, a)
    xs = np.linspace(0, math.pi, n, endpoint=False)
    good = np.ones(n, bool)
    for s in S:
        lhs = np.abs(np.sin(xs - s*math.pi/4))
        rhs = abs(math.sin(th - s*math.pi/4))
        good &= lhs >= rhs - 1e-12
    return xs[good], th

bad = 0; cnt = 0
for p in split_primes(5, 20000):
    a, b = gauss_prime(p)
    S = [2, 3] if b < a else [1, 2, 3]
    xs, th = sols(a, b, S, n=40000)
    cands = [th] if b < a else [th, math.pi/2 - th]
    h = math.pi/40000
    # every grid solution must be within 2 grid steps of a candidate
    for x in xs:
        if min(abs(x - c) for c in cands) > 2.5*h:
            bad += 1; break
    cnt += 1
print(f"(b) grid check (40000 pts mod pi) on {cnt} split primes <= 2e4: grid solutions off the candidate set: {bad}")

# (c) p-part condition (i): no p_l divides a_j, a_j+b_j, and (if b_j>a_j) a_j-b_j
def fail_i(ps):
    P = set(ps); f = 0
    for p in ps:
        a, b = gauss_prime(p)
        nums = [a, a+b] + ([abs(a-b)] if b > a else [])
        if any(any(l in P for l in factorint(x)) for x in nums if x > 1):
            f += 1
    return f

random.seed(14)
print("(c) columns failing condition (i) [p_l | a_j(a_j+b_j)(a_j-b_j)], random prime sets:")
for (lo, hi, r) in [(10**2, 10**3, 50), (10**3, 10**4, 100), (10**4, 10**5, 200), (10**5, 10**6, 200), (10**2, 10**6, 200)]:
    pool = split_primes(lo, hi)
    tot = 0; T = 20
    for _ in range(T):
        ps = random.sample(pool, r)
        tot += fail_i(ps)
    print(f"   primes in [{lo},{hi}], r={r}: mean failing columns {tot/T:.2f} of {r}")
