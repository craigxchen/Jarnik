# Independent checks for the "meaning" block (own code, exact integers).
import random, math
from math import gcd, isqrt
from fractions import Fraction
from sympy import factorint
random.seed(31)

# (1) Square-block identity of VIII.2 on random two-class CRT 8-tuples (with coalescing residues).
def crt(rs, ms):
    x, m = 0, 1
    for r, mm in zip(rs, ms):
        # solve x = x0 mod m, x = r mod mm
        t = ((r - x) * pow(m, -1, mm)) % mm
        x += m * t; m *= mm
    return x % m, m
primes = [3,5,7,11,13,17,19,23,29,31,37,41,43]
fails = 0; checks = 0
for trial in range(3000):
    ps = random.sample(primes, random.randint(1,5))
    mods = [p**random.randint(1,3) for p in ps]
    q = 1
    for m in mods: q *= m
    classes = []
    for p, m in zip(ps, mods):
        r1 = random.randrange(m)
        e = round(math.log(m, p)); s = random.randint(0, e-1)
        r2 = (r1 + p**s * random.randrange(1, p)) % m  # differ exactly at level s
        if random.random() < 0.2: r2 = r1
        classes.append((r1, r2))
    # 8 distinct marks: random choice of class per prime + random lift
    marks = set()
    while len(marks) < 8:
        rs = [c[random.randrange(2)] for c in classes]
        x, _ = crt(rs, mods)
        marks.add(x + q * random.randint(-6, 6))
    X = list(marks); random.shuffle(X)
    xs = X[0]
    g = q
    for x in X: g = gcd(g, x - xs)
    Q = q // g
    Y = [(x - xs) // g for x in X]
    if Q == 1: continue
    fQ = factorint(Q)
    y0 = Y[0]
    d = [gcd(Q, y - y0) for y in Y]
    P = [(Q // d[i], (Y[i] - y0) // d[i]) for i in range(8)]
    assert P[0] == (1, 0)
    for i in range(8): assert gcd(P[i][0], P[i][1]) == 1
    # blocks: for each p^f || Q, T_p = slots opposite anchor mod p^f
    blocks = {}
    for p, f in fQ.items():
        pf = p**f
        T = frozenset(i for i in range(1, 8) if (Y[i] - y0) % pf)
        # check two classes with unit difference
        cls = set(y % pf for y in Y); assert len(cls) <= 2
        if len(cls) == 2:
            a, b = cls; assert (a - b) % p != 0
        if T: blocks[T] = blocks.get(T, 1) * p**(2*f)
    for i in range(8):
        for j in range(i+1, 8):
            det = P[i][0]*P[j][1] - P[i][1]*P[j][0]
            prod = 1
            for T, n in blocks.items():
                if (i == 0 or i in T) and j in T: prod *= n  # anchor is never in T; T contains i,j
            # convention: T subset of [7]; minor (0,j) involves blocks containing j and 0 -> none
            prod = 1
            for T, n in blocks.items():
                if i in T and j in T: prod *= n
            h = gcd(Q, Y[j] - Y[i])
            t = (Y[j] - Y[i]) // h
            checks += 1
            if det != t * prod: fails += 1
            # cofactor size identity: log|t| = D/2 - W/4 + log(|yj-yi|/sqrt Q)
            D = Fraction(Q*Q, h*h)  # e^{D}
            # |t|^2 = e^D * (yj-yi)^2 / Q^2 * ... check: t^2 = (yj-yi)^2/h^2 ; e^{D/2 - W/4}=(Q/h)/sqrt(Q)
            assert Fraction(t*t) == D * Fraction((Y[j]-Y[i])**2, Q*Q)
print("(1) square-block identity: pair checks", checks, "failures", fails)

# (2) Liouville exponents: Cilleruelo-Cordoba arcs vs CG Theorem 10 (two residues)
ok = True
for k in range(1, 60):
    cc = Fraction(1,2) - Fraction(1, 4*(k//2)+2)      # k+1 points force arc length >= R^cc
    kk = k + 1; l = kk if kk % 2 else kk - 1           # k+1 roots force span > b^(1/2-1/(2l))
    cg = Fraction(1,2) - Fraction(1, 2*l)
    if cc != cg: ok = False; print("differ", k, cc, cg)
print("(2) CC exponent == CG Thm10 exponent for k=1..59:", ok)
