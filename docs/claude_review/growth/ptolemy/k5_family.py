r"""Exact verification of the explicit k=5 bounded-residue balanced Pluecker family (Theorem C1).

Rows P_0=(1,0), P_i=(X_i(s),1), i=1..4, with
    X_i(s) = a_i (s - tau_0) prod_{l != i} (s - tau_l)          (degree 4)
Blocks (T subset [4], |T|>=2):
    n_[4]        = s - tau_0                         (top block)
    n_{[4]\{l}}  = s - tau_l                         (triples)
    n_ij         = sign_ij * ((a_i-a_j) s - (a_i tau_j - a_j tau_i))   (pairs)
Claims (checked below exactly, as polynomial identities in s and at many integer s):
  (1) det(P_0,P_i) = 1, det(P_i,P_j) = t_ij prod_{T contains i,j} n_T with t_ij = +-1;
  (2) every one of the C(5,4)=5 quadruple relations holds in the matching form
        n_A a + n_B b = n_C c     with |a|,|b|,|c| = 1  (signs from orientation);
  (3) all 11 blocks are pairwise coprime for s in an explicit residue class, and
      log n_T = log s + O(1) for every T (balanced profile, w = log s);
  (4) the cross-ratio valuations at every prime are exactly the Boolean ones (p-adic tree
      edges = blocks), i.e. the M_{0,5}-point meets boundary divisor D_T exactly at the primes of n_T.
"""
import itertools
from fractions import Fraction
from math import gcd, log

A = [None, 1, 2, 3, 4]           # a_1..a_4
TAU = [4, 6, 7, 0, 2]            # tau_0 (top), tau_1..tau_4 (triple [4]\{l})

def polymul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            r[i + j] += x * y
    return r

def polysub(p, q):
    n = max(len(p), len(q)); p = [0] * (n - len(p)) + p; q = [0] * (n - len(q)) + q
    r = [x - y for x, y in zip(p, q)]
    while len(r) > 1 and r[0] == 0:
        r.pop(0)
    return r

def polyval(p, s):
    v = 0
    for c in p:
        v = v * s + c
    return v

def lin(alpha, beta):          # alpha*s - beta, coefficients highest first
    return [alpha, -beta]

ROWS = [1, 2, 3, 4]
def X(i):
    p = [A[i]]
    p = polymul(p, lin(1, TAU[0]))
    for l in ROWS:
        if l != i:
            p = polymul(p, lin(1, TAU[l]))
    return p

blocks = {}
blocks[frozenset(ROWS)] = lin(1, TAU[0])
for l in ROWS:
    blocks[frozenset(x for x in ROWS if x != l)] = lin(1, TAU[l])
for i, j in itertools.combinations(ROWS, 2):
    blocks[frozenset((i, j))] = lin(A[i] - A[j], A[i] * TAU[j] - A[j] * TAU[i])

def det_poly(i, j):
    """det(P_i,P_j) as polynomial; P_0=(1,0), P_i=(X_i,1)"""
    if i == 0:
        return [1]
    if j == 0:
        return [-1]
    return polysub(X(i), X(j))

def G(i, j):
    p = [1]
    for T, b in blocks.items():
        if i in T and j in T:
            p = polymul(p, b)
    return p

# (1) exact factorization with residual +-1 (as polynomials)
t = {}
for i, j in itertools.combinations([0] + ROWS, 2):
    d = det_poly(i, j)
    if i == 0:
        assert d == [1]; t[i, j] = 1; continue
    g = G(i, j)
    if d == g:
        t[i, j] = 1
    elif d == [-c for c in g]:
        t[i, j] = -1
    else:
        raise AssertionError((i, j, d, g))
print("(1) det(P_i,P_j) = t_ij * prod_{T>=ij} n_T with t_ij in {+1,-1}:", t)

# roots distinct
roots = {T: Fraction(b[1], -b[0]) * -1 if False else Fraction(-b[1], b[0]) for T, b in blocks.items()}
assert len(set(roots.values())) == len(roots), roots
print("    the 11 block roots are distinct:", sorted(roots.values()))

# (2) all quadruple relations in matching form, as polynomial identities
def p_(i, j):
    return det_poly(i, j) if i < j else [-c for c in det_poly(j, i)]
cnt = 0
for Q in itertools.combinations([0] + ROWS, 4):
    a, b, c, d = Q
    # Pluecker: p_ac p_bd = p_ab p_cd + p_ad p_bc
    lhs = polymul(p_(a, c), p_(b, d))
    rhs = polysub(polymul(p_(a, b), p_(c, d)), [-x for x in polymul(p_(a, d), p_(b, c))])
    assert polysub(lhs, rhs) == [0]
    # matching coefficients: blocks splitting Q 2|2 (side not containing anchor 0)
    def nprod(pairs):
        p = [1]
        for T, bl in blocks.items():
            if frozenset(Q) & T in [frozenset(x) for x in pairs] and 0 not in T:
                p = polymul(p, bl)
        return p
    nA = nprod([(a, b), (c, d)]); nB = nprod([(a, d), (b, c)]); nC = nprod([(a, c), (b, d)])
    # residual products
    tt = lambda u, v: t[min(u, v), max(u, v)] * (1 if u < v else -1)
    ra, rb, rc = tt(a, b) * tt(c, d), tt(a, d) * tt(b, c), tt(a, c) * tt(b, d)
    # remove the common factor: p_xy = t_xy * G_xy ; G_ab G_cd = common * nA etc.
    common = None
    lhs2 = polymul([rc], nC); rhs2 = polysub(polymul([ra], nA), [-x for x in polymul([rb], nB)])
    assert polysub(lhs2, rhs2) == [0], (Q, nA, nB, nC)
    cnt += 1
    print(f"(2) quadruple {Q}: n_C*({rc}) = n_A*({ra}) + n_B*({rb}),  deg n_A,n_B,n_C =",
          len(nA) - 1, len(nB) - 1, len(nC) - 1)
assert cnt == 5

# (3) coprimality class and balance.  Primes q>=5 dividing a resultant: choose s mod q with no
# block divisible by q.  Primes 2,3: choose s mod q^E so that each block has a CONSTANT q-valuation
# (< E) on the class; these bounded q-parts c_T are moved into the residues.
def factor(n):
    n = abs(n); f = set(); p = 2
    while p * p <= n:
        while n % p == 0:
            f.add(p); n //= p
        p += 1
    if n > 1:
        f.add(n)
    return f
def vq(n, q):
    n = abs(n); v = 0
    while n % q == 0:
        n //= q; v += 1
    return v
bl = list(blocks.items())
bad_primes = {2, 3}
for (T1, b1), (T2, b2) in itertools.combinations(bl, 2):
    res = b1[0] * b2[1] - b2[0] * b1[1]
    assert res != 0
    bad_primes |= factor(res)
cls = {}
E = 4
for q in sorted(bad_primes):
    if q >= 5:
        ok = [r for r in range(q) if all(polyval(b, r) % q != 0 for _, b in bl)]
        assert ok, ("no admissible class mod", q)
        cls[q] = ok[0]
    else:
        ok = [r for r in range(q ** E) if all(polyval(b, r) % (q ** E) != 0 for _, b in bl)]
        assert ok
        cls[q ** E] = ok[0]
Q = 1; s0 = 0
for q, r in cls.items():
    while s0 % q != r:
        s0 += Q
    Q *= q
print(f"(3) residue class s = {s0} mod {Q} (moduli {sorted(cls)})")
def reduced(T, s):
    v = abs(polyval(blocks[T], s))
    c = 2 ** vq(v, 2) * 3 ** vq(v, 3)
    return v // c, c
checked = 0; cvals = None
for k in range(1, 400):
    s = s0 + Q * k
    red = {T: reduced(T, s) for T in blocks}
    cv = {T: c for T, (v, c) in red.items()}
    if cvals is None:
        cvals = cv
    assert cv == cvals                      # the {2,3}-parts are constant on the class
    for (T1, (v1, _)), (T2, (v2, _)) in itertools.combinations(red.items(), 2):
        assert gcd(v1, v2) == 1
    dev = max(abs(log(v) - log(s)) for v, _ in red.values())
    assert dev < log(20)
    checked += 1
tmax = max(abs(t[i, j]) * __import__('math').prod(cvals[T] for T in blocks if i in T and j in T)
           for i, j in itertools.combinations(ROWS, 2))
print(f"    constant {{2,3}}-parts of the blocks on the class: {sorted(set(cvals.values()))};")
print(f"    with them moved into the residues, |t'_ij| <= {tmax}, |Y_i| = 1;")
print(f"    reduced blocks pairwise coprime and |log n'_T - log s| < log 20, verified at {checked} values of s")

# (4) p-adic tree: cross-ratio valuations equal Boolean values at every prime, for sample s
def vp(n, p):
    n = abs(n); v = 0
    while n % p == 0:
        n //= p; v += 1
    return v
for k in (3, 17, 101):
    s = s0 + Q * k
    P = {0: (1, 0)}
    for i in ROWS:
        P[i] = (polyval(X(i), s), 1)
    det = lambda u, v: P[u][0] * P[v][1] - P[v][0] * P[u][1]
    primes = set()
    for T, b in bl:
        primes |= factor(polyval(b, s))
    primes -= {2, 3}
    for p in primes:
        for (a, b2, c, d) in itertools.permutations([0] + ROWS, 4):
            if not (a < b2 and c < d and a < c):
                continue
            # valuation of cross-ratio [ab|cd] = det(a,b)det(c,d)/(det(a,c)det(b,d))
            v = vp(det(a, b2), p) + vp(det(c, d), p) - vp(det(a, c), p) - vp(det(b2, d), p)
            boolv = 0
            for T, bb in bl:
                vt = vp(polyval(bb, s), p)
                if vt == 0:
                    continue
                ins = lambda x: x in T
                boolv += vt * (int(ins(a) and ins(b2)) + int(ins(c) and ins(d)) - int(ins(a) and ins(c)) - int(ins(b2) and ins(d)))
            assert v == boolv
print("(4) all cross-ratio valuations at all block primes agree with the Boolean tree (3 sample s)")
print("ALL CHECKS PASSED")
