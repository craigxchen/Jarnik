# Judge check (i): IX.0 main theorem. Own code, exact integers.
# (a) union-bound exponent recomputed with exact rationals / floats for logs
# (b) conditional law: exact enumeration on a small pool that the hidden k-tuple, given exposed slots,
#     is uniform over injective k-tuples of unused primes, and the bound 1/binom(K-s+k,k) <= (s/K)^k
# (c) twin: for random squarefree q, all admissible triples (0,a,b) in [0,L]: identity (IX.0a),
#     exactly one block per non-constant prime, heavy block >= (W/4 - log C/2)/log maxp primes
# (d) Gaussian: random points of norm q (split primes), same unit: z_j - z_i = 2 i eta H t, triangle in Z[i]
import math, random, itertools
from fractions import Fraction
from math import comb, log

random.seed(20261120)
out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

# (a)
for labels, name in ((4, "twin"), (7, "gauss")):
    for J in (13, 18.6, 20, 40, 60):
        per = math.log(labels)/J + 9/200 - 1/5   # per sJ, b-1 = 100 log 2P / J
        say(f"(a) {name} J={J}: exponent per sJ = {per:.4f}  (need <= -0.05: {per <= -0.05})")
# bounded-exponent twin: labels 1+3E, k >= s/(5E), coefficient 9E s log2P/(2(b-1)), b-1 >= 100E^2 log2P/J
for E in (1, 2, 3, 5):
    J = 40*E*math.log(1+3*E)
    per = math.log(1+3*E)/J + 9*E/(2*100*E*E) - 1/(5*E)  # per sJ
    say(f"(a) twin powers E={E}, J={J:.1f}: exponent per sJ = {per:.4f} = {per*E:.4f}/E")

# (b) exact conditional law on a toy pool
def toy(K, s, k):
    pool = list(range(K))
    # fix labelling: first s-k slots exposed, last k slots hidden (any fixed labelling is equivalent by symmetry
    # because the order is uniform); enumerate all ordered injective s-tuples
    from collections import Counter
    cond = {}
    for seq in itertools.permutations(pool, s):
        e = seq[:s-k]; h = frozenset(seq[s-k:])
        cond.setdefault(e, Counter())[h] += 1
    worst = 0
    for e, c in cond.items():
        tot = sum(c.values())
        assert len(c) == comb(K-(s-k), k)
        for h, v in c.items():
            worst = max(worst, Fraction(v, tot))
    bound = Fraction(1, comb(K-s+k, k))
    return worst, bound, Fraction(s, K)**k
for (K, s, k) in ((7, 3, 2), (8, 4, 2), (8, 4, 3), (9, 3, 1)):
    w, b, r = toy(K, s, k)
    say(f"(b) K={K} s={s} k={k}: max P(hidden set = given set | exposed) = {w}, 1/binom = {b}, (s/K)^k = {r}; ok={w == b and b <= r}")
# inequality binom(K-s+k,k)^-1 <= (s/K)^k over a grid
bad = 0
for K in range(2, 120):
    for s in range(1, K+1):
        for k in range(0, s+1):
            if Fraction(1, comb(K-s+k, k)) > Fraction(s, K)**k: bad += 1
say("(b) grid K<120: violations of 1/binom(K-s+k,k) <= (s/K)^k:", bad)

# (c) twin triples
def primes_in(a, b):
    ps = []
    for n in range(max(2, a), b+1):
        if all(n % p for p in range(2, int(n**0.5)+1)): ps.append(n)
    return ps
def twin_check(trials=40):
    nt = 0; nfail = 0; minslack = None
    for _ in range(trials):
        s = random.randint(3, 5)
        ps = random.sample(primes_in(30, 90), s); q = math.prod(ps)
        C = random.choice([1, 2, 3]); L = int(C*math.isqrt(q))
        W = log(q)
        # admissible triples (0,a,b): q | a*b*(b-a)  (squarefree q: at most two classes per p)
        for a in range(1, L+1):
            g = math.gcd(a, q)
            r = q//g
            rp = [p for p in ps if r % p == 0]
            bs = set()
            for mask in range(1 << len(rp)):
                r1 = math.prod(rp[i] for i in range(len(rp)) if mask >> i & 1); r2 = r//r1
                # b = 0 mod r1, b = a mod r2
                b0 = (a * r1 * pow(r1, -1, r2)) % r if r2 > 1 else 0
                b0 %= r
                bb = b0
                while bb <= a: bb += r
                while bb <= L: bs.add(bb); bb += r
            for b in sorted(bs):
                y = (0, a, b); nt += 1
                h = lambda i, j: math.prod(p for p in ps if (y[j]-y[i]) % p == 0)
                G = math.prod(p for p in ps if a % p == 0 and b % p == 0)
                A12 = h(0,1)//G; A23 = h(1,2)//G; A13 = h(0,2)//G
                t12 = a//h(0,1); t23 = (b-a)//h(1,2); t13 = b//h(0,2)
                ok = (t12*A12 + t23*A23 == t13*A13) and G*A12*A23*A13 == q
                ok &= all(math.gcd(x, z) == 1 for x, z in ((G,A12),(G,A23),(G,A13),(A12,A23),(A12,A13),(A23,A13)))
                # heavy block on this triple: separating 1,2 = A13*A23
                big = max(A13, A23)
                cnt = sum(1 for p in ps if big % p == 0)
                need = (W/4 - log(C)/2)/log(max(ps))
                if cnt < need - 1e-9: ok = False
                if not ok: nfail += 1
    return nt, nfail
nt, nf = twin_check()
say(f"(c) twin: {nt} admissible triples checked, failures (identity, coprime factorisation, heavy-block count) = {nf}")

# (d) Gaussian
def gmul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gconj(x): return (x[0], -x[1])
def gpow(x, e):
    r = (1, 0)
    for _ in range(e): r = gmul(r, x)
    return r
def gen(p):
    for a in range(1, int(p**0.5)+1):
        b2 = p - a*a; b = math.isqrt(b2)
        if b*b == b2 and a > b: return (a, b)
    raise ValueError
def gdivexact(x, y):
    n = y[0]**2 + y[1]**2; num = gmul(x, gconj(y))
    assert num[0] % n == 0 and num[1] % n == 0, "not divisible"
    return (num[0]//n, num[1]//n)
split = [p for p in primes_in(5, 400) if p % 4 == 1]
nf = 0; ntr = 0
for _ in range(300):
    s = random.randint(3, 7); ps = random.sample(split, s); pis = {p: gen(p) for p in ps}
    eta = random.choice([(1,0),(0,1),(-1,0),(0,-1)])
    marks = []
    for _ in range(3):
        g = {p: random.randint(0, 1) for p in ps}
        marks.append(g)
    if marks[0] == marks[1] or marks[1] == marks[2] or marks[0] == marks[2]: continue
    def z(g):
        r = eta
        for p in ps: r = gmul(r, pis[p] if g[p] == 0 else gconj(pis[p]))
        return r
    Z = [z(g) for g in marks]
    def H(i, j):
        r = (1, 0)
        for p in ps:
            if marks[i][p] == marks[j][p]: r = gmul(r, pis[p] if marks[i][p] == 0 else gconj(pis[p]))
        return r
    ok = True; T = {}
    for (i, j) in ((0,1),(1,2),(0,2)):
        d = (Z[j][0]-Z[i][0], Z[j][1]-Z[i][1])
        tq = gdivexact(d, gmul((0, 2), gmul(eta, H(i, j))))
        if tq[1] != 0 or tq[0] == 0: ok = False
        T[(i, j)] = tq[0]
    # triangle after dividing by common G
    Gc = (1, 0)
    for p in ps:
        if marks[0][p] == marks[1][p] == marks[2][p]: Gc = gmul(Gc, pis[p] if marks[0][p] == 0 else gconj(pis[p]))
    A = {k: gdivexact(H(*k), Gc) for k in T}
    lhs = tuple(T[(0,1)]*A[(0,1)][c] + T[(1,2)]*A[(1,2)][c] for c in (0,1))
    rhs = tuple(T[(0,2)]*A[(0,2)][c] for c in (0,1))
    if lhs != rhs: ok = False
    ntr += 1; nf += (not ok)
say(f"(d) gaussian squarefree: {ntr} triples, failures = {nf}")
open(__file__.replace('.py', '.out'), 'w').write("\n".join(out) + "\n")
