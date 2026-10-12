"""Checks for the height/valuation lemmas of lw.md, Section 2 (exact integer arithmetic).

Part A: for random alpha in Q(i)^* (non-torsion), the primitive minimal polynomial over Z, its naive
        height H, and P^+(alpha), P^-(alpha) (sums of ord_v log N(v) over Gaussian primes with positive /
        negative valuation).  Checks: M(alpha) = exp(max(P^+, P^-)) for quadratic alpha (= 2h),
        H(alpha) = exp(max(P^+,P^-)/2) for rational alpha, and M <= H <= 2 M (quadratic; Lemma 2.1(iii)).
Part B: brute force over representations  T = prod_k alpha_k^{c_k} * unit  with K <= 2 numbers from a
        finite pool (all products pi_1^x1 conj(pi_1)^y1 pi_2^x2 conj(pi_2)^y2, exponents in [-2,2]) and
        |c_k| <= 4, for targets T = u_1^{v1} u_2^{v2} (u_j = pi_j/conj pi_j): the minimum of
        sum_k log H(alpha_k) is compared with the support weight w(supp v) of Lemma 2.3.
Usage: python3 check_valuation.py
"""
import math, random, itertools
from fractions import Fraction

random.seed(3)
out = []
def log(s):
    print(s); out.append(s)

# Gaussian integers as tuples (a,b) = a+bi
def gmul(x, y):
    return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def gconj(x):
    return (x[0], -x[1])
def gnorm(x):
    return x[0]*x[0] + x[1]*x[1]
def gpow(x, e):
    r = (1, 0)
    for _ in range(e):
        r = gmul(r, x)
    return r

def gdivides(p, x):
    # does Gaussian integer p divide x?
    n = gnorm(p)
    q = gmul(x, gconj(p))
    return q[0] % n == 0 and q[1] % n == 0
def gdiv(x, p):
    n = gnorm(p)
    q = gmul(x, gconj(p))
    return (q[0] // n, q[1] // n)

def gaussian_primes_upto(B):
    ps = []
    for p in range(2, B + 1):
        if all(p % d for d in range(2, int(p**0.5) + 1)):
            if p == 2:
                ps.append(((1, 1), 2))
            elif p % 4 == 1:
                for a in range(1, int(p**0.5) + 1):
                    b2 = p - a*a
                    b = int(round(b2**0.5))
                    if b*b == b2:
                        ps.append(((a, b), p)); ps.append(((a, -b), p)); break
            else:
                ps.append(((p, 0), p*p))
    return ps
GP = gaussian_primes_upto(200)

def valuations(x):
    """dict prime-index -> ord, for a Gaussian integer with all prime factors in GP (norm <= 200^2)."""
    v = {}
    for i, (p, n) in enumerate(GP):
        while gnorm(x) > 1 and gdivides(p, x):
            x = gdiv(x, p); v[i] = v.get(i, 0) + 1
    assert gnorm(x) == 1, x
    return v

def minpoly(num, den):
    """primitive minimal polynomial over Z of alpha = num/den (Gaussian integers)."""
    # alpha = num*conj(den)/N(den) = (A + B i)/D
    t = gmul(num, gconj(den)); D = gnorm(den)
    A, B = t
    g = math.gcd(math.gcd(A, B), D)
    A //= g; B //= g; D //= g
    if B == 0:
        g2 = math.gcd(A, D)
        return [D // g2, -A // g2]  # D x - A
    # (D x - A)^2 + B^2 = D^2 x^2 - 2AD x + A^2 + B^2
    c = [D*D, -2*A*D, A*A + B*B]
    g3 = math.gcd(math.gcd(c[0], c[1]), c[2])
    return [x // g3 for x in c]

def mahler(poly):
    if len(poly) == 2:
        return max(abs(poly[0]), abs(poly[1]))
    a, b, c = poly
    # complex conjugate roots of modulus sqrt(c/a)
    r2 = Fraction(c, a)
    return abs(a) * max(1, r2)  # |a| * max(1,|r|)^2

# Part A
nA = 0
for trial in range(3000):
    # random numerator / denominator from small Gaussian primes
    def rand_g():
        x = (1, 0)
        for _ in range(random.randint(0, 3)):
            p, n = random.choice(GP[:25])
            x = gmul(x, p)
        return x
    num, den = rand_g(), rand_g()
    if num == (0, 0):
        continue
    poly = minpoly(num, den)
    vn, vd = valuations(num), valuations(den)
    ordv = {}
    for k, e in vn.items(): ordv[k] = ordv.get(k, 0) + e
    for k, e in vd.items(): ordv[k] = ordv.get(k, 0) - e
    Pp = sum(e * math.log(GP[k][1]) for k, e in ordv.items() if e > 0)
    Pm = sum(-e * math.log(GP[k][1]) for k, e in ordv.items() if e < 0)
    H = max(abs(x) for x in poly)
    Mh = mahler(poly)
    if len(poly) == 3:
        if Pp == 0 and Pm == 0:
            continue  # torsion (unit)
        assert abs(math.log(Mh) - max(Pp, Pm)) < 1e-9, (num, den, poly)
        assert Mh - 1e-9 <= H <= 2 * Mh + 1e-9  # H >= M exactly for quadratic alpha in Q(i)
    else:
        if Pp == 0 and Pm == 0:
            continue
        assert abs(math.log(H) - max(Pp, Pm) / 2) < 1e-9, (num, den, poly)
    nA += 1
log(f"Part A: {nA} random alpha in Q(i): 2h(alpha) = max(P+,P-) (quadratic), log H = max(P+,P-)/2 (rational), "
    f"M <= H <= 2M (quadratic): all PASS")

# Part B
pi1, p1 = (2, 1), 5
pi2, p2 = (3, 2), 13
def elt(x1, y1, x2, y2):
    """alpha = pi1^x1 conj(pi1)^y1 pi2^x2 conj(pi2)^y2 as (num, den) Gaussian integers."""
    num, den = (1, 0), (1, 0)
    for base, e in ((pi1, x1), (gconj(pi1), y1), (pi2, x2), (gconj(pi2), y2)):
        if e >= 0: num = gmul(num, gpow(base, e))
        else: den = gmul(den, gpow(base, -e))
    return num, den
pool = []
for ex in itertools.product(range(-2, 3), repeat=4):
    if ex == (0, 0, 0, 0):
        continue
    num, den = elt(*ex)
    poly = minpoly(num, den)
    pool.append((ex, math.log(max(abs(c) for c in poly)), len(poly) == 3))
w1, w2 = math.log(p1), math.log(p2)
for v in ((1, 0), (1, 1), (1, -1), (2, -1), (2, 1), (3, -2)):
    # target exponent vector on (pi1, conj pi1, pi2, conj pi2): u_j^{v_j} = pi_j^{v_j} conj(pi_j)^{-v_j}
    tgt = (v[0], -v[0], v[1], -v[1])
    best = (1e9, None)
    for K in (1, 2):
        for combo in itertools.combinations(range(len(pool)), K):
            for cs in itertools.product([c for c in range(-4, 5) if c != 0], repeat=K):
                tot = [sum(cs[k] * pool[combo[k]][0][i] for k in range(K)) for i in range(4)]
                if tuple(tot) == tgt:
                    cost = sum(pool[combo[k]][1] for k in range(K))
                    if cost < best[0]:
                        best = (cost, [(pool[combo[k]][0], cs[k]) for k in range(K)])
    wS = (w1 if v[0] else 0) + (w2 if v[1] else 0)
    nq = 2
    log(f"Part B target v={v}: min sum log H over representations (K<=2, |c|<=4) = {best[0]:.4f}; "
        f"w(supp v) = {wS:.4f}; Lemma 2.3 lower bound w(supp) = {wS:.4f}; best rep {best[1]}")
    assert best[0] >= wS - 1e-9
open('check_valuation_out.txt', 'w').write('\n'.join(out) + '\n')
