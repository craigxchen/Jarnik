"""Top degree of Syl(16, L) (T = 1 Sylvester box, n = 16) via the Singer cycle and k -> -k.

Coordinates: nonzero points alpha^i of F_16 (alpha root of x^4 + x + 1), i = 0..14; base point 0.
Columns: S_c = {i : Tr(c alpha^i) = 1}, c != 0: the 15 cyclic shifts of D = {i : Tr(alpha^i) = 1}.
Group Z_15 x Z_2 (shift, negation) acts on A; for each character (s, eps) the isotypic component of the function
space has a basis indexed by orbit representatives with compatible stabiliser, and is spanned by the projections of
monomial-orbit representatives of parity eps.  We test whether the monomials of degree <= j (j = argv[2]) span
every component: if all components have full rank mod p (prime p = 1 mod 15), then rho_Q <= rho_p <= j.
Columns are a random sample of monomial orbits of degree <= j and right parity (sampling can only lower the rank,
so 'full rank' is a valid certificate). In this version ALL orbits of degree <= j are used, so the answer is exact mod p.
"""
import sys, time, itertools, random, os
import numpy as np
import flint
from sympy import isprime

q = 15
L = int(sys.argv[1]); jdeg = int(sys.argv[2]); extra = float(sys.argv[3]) if len(sys.argv) > 3 else 1.15
t0 = time.time()

# F_16 trace pattern
def gf_mul_alpha(v):  # v = 4-bit int, multiply by alpha mod x^4+x+1
    v <<= 1
    if v & 16:
        v ^= 0b10011
    return v
powers = [1]
for i in range(14):
    powers.append(gf_mul_alpha(powers[-1]))
def trace(v):
    t = 0; x = v
    for _ in range(4):
        t ^= x
        # square in GF(16)
        # compute x^2 by multiplying
        y = 0; a = x; b = x
        for bit in range(4):
            if (b >> bit) & 1:
                aa = a
                for _ in range(bit):
                    aa = gf_mul_alpha(aa)
                y ^= aa
        x = y
    # t is an element of F_2 embedded (0 or 1)
    return t & 1
D = [i for i in range(q) if trace(powers[i]) == 1]
assert len(D) == 8, D
cols = [sorted((i + t) % q for i in D) for t in range(q)]

# enumerate A with numpy
C = np.zeros((q, q), dtype=np.int64)
for p, S in enumerate(cols):
    C[p, S] = 1
first = np.array(list(itertools.product((-1, 0, 1), repeat=8)), dtype=np.int64)
rest = np.array(list(itertools.product((-1, 0, 1), repeat=q - 8)), dtype=np.int64)
vf = first @ C[:, :8].T
vr = rest @ C[:, 8:].T
pts = []
for i in range(len(rest)):
    ok = np.all(np.abs(vf + vr[i]) <= L, axis=1)
    sel = first[ok]
    if len(sel):
        pts.append(np.hstack([sel, np.repeat(rest[i:i + 1], len(sel), axis=0)]))
P = np.vstack(pts)
N = len(P)
print('Syl16 L=%d |A|=%d D=%s (%.0fs)' % (L, N, D, time.time() - t0), flush=True)

# orbit representatives under shift x negation
pw3 = 3 ** np.arange(q, dtype=np.int64)
def code(M):
    return ((M + 1) * pw3).sum(axis=1)
codes = []
for g in range(q):
    Sg = np.roll(P, g, axis=1)
    codes.append(code(Sg)); codes.append(code(-Sg))
codes = np.array(codes)  # 30 x N
canon = codes.min(axis=0)
own = code(P)
is_rep = canon == own
R = P[is_rep]
# stabiliser info: for each rep, list of (g, e) fixing it
Rcodes = own[is_rep]
stab = []
for idx in range(len(R)):
    st = []
    for g in range(q):
        for e, sgn in ((0, 1), (1, -1)):
            if codes[2 * g + e, np.where(is_rep)[0][idx]] == Rcodes[idx]:
                st.append((g, e))
    stab.append(st)
print('orbit reps %d (%.0fs)' % (len(R), time.time() - t0), flush=True)

# prime p = 1 mod 15, < 2^59
p = (1 << 58) - ((1 << 58) % 15) + 1
while not isprime(p):
    p -= 15
rng = random.Random(7)
h = 2
while True:
    om = pow(h, (p - 1) // 15, p)
    if pow(om, 5, p) != 1 and pow(om, 3, p) != 1:
        break
    h += 1

# monomial orbit reps of degree <= jdeg (exponents in {0,1,2}), by parity; enumerate canonical under shift
def mon_orbits(maxdeg):
    reps = {0: [], 1: []}
    for a in itertools.product((0, 1, 2), repeat=q):
        dg = sum(a)
        if dg > maxdeg:
            continue
        rots = [a[i:] + a[:i] for i in range(q)]
        if a != min(rots):
            continue
        reps[dg % 2].append(a)
    return reps
MO = mon_orbits(jdeg)
print('monomial orbits deg<=%d: even %d odd %d (%.0fs)' % (jdeg, len(MO[0]), len(MO[1]), time.time() - t0), flush=True)

pw2 = (1 << np.arange(q, dtype=np.int64))
def masks_pts(M):
    nz = ((M != 0).astype(np.int64) * pw2).sum(axis=1)
    ng = ((M < 0).astype(np.int64) * pw2).sum(axis=1)
    return nz, ng
PARTAB = np.array([bin(i).count('1') & 1 for i in range(1 << q)], dtype=np.int64)
def parity_pop(x):
    return PARTAB[x]

allfull = True
HIGH = int(os.environ.get('HIGH', '800'))
ONLY = [tuple(map(int, x.split(','))) for x in os.environ.get('ONLY', '').split(';') if x]
for eps in (0, 1):
    monos = MO[eps]
    for s in range(q):
        if ONLY and (eps, s) not in ONLY:
            continue
        # compatible reps: character chi(g,e) = om^{s g} * (-1)^{eps e}; need chi = 1 on stabiliser
        rows = [i for i in range(len(R)) if all((pow(om, (s * g) % q, p) * (1 if (eps * e) % 2 == 0 else p - 1)) % p == 1 for g, e in stab[i])]
        nr = len(rows)
        if nr == 0:
            continue
        # column plan (sampling can only lower the rank, so FULL rank is a certificate that degree <= jdeg spans)
        if eps == 0:
            plan = [((0, 8), None), ((10, jdeg), HIGH)]
        else:
            plan = [((0, 7), None), ((9, 9), 4500), ((11, jdeg), HIGH)]
        sample = []
        for (lo_, hi_), cnt in plan:
            pool = [a for a in monos if lo_ <= sum(a) <= hi_]
            sample += pool if (cnt is None or cnt >= len(pool)) else rng.sample(pool, cnt)
        A1 = np.array([[1 if x == 1 else 0 for x in a] for a in sample], dtype=np.int64)
        A2 = np.array([[1 if x == 2 else 0 for x in a] for a in sample], dtype=np.int64)
        a1 = (A1 * pw2).sum(axis=1); asupp = ((A1 + A2) * pw2).sum(axis=1)
        acc = np.zeros((nr, len(sample)), dtype=np.int64)
        Rr = R[rows]
        for g in range(q):
            # value of mu at T^g k : np.roll(k, g) means (T^g k)_j = k_{j-g}
            Sg = np.roll(Rr, g, axis=1)
            nz, ng = masks_pts(Sg)
            zero = (asupp[None, :] & ~nz[:, None]) != 0
            par = parity_pop(a1[None, :] & ng[:, None])
            val = np.where(zero, 0, 1 - 2 * par).astype(np.int64)
            coef = pow(om, (-s * g) % q, p)
            acc = (acc + val * coef) % p
        M = flint.nmod_mat(nr, len(sample), acc.ravel().tolist(), p)
        rk = M.rank()
        full = (rk == nr)
        allfull &= full
        print('eps=%d s=%d rows=%d cols=%d rank=%d %s (%.0fs)' % (eps, s, nr, len(sample), rk, 'FULL' if full else 'deficient', time.time() - t0), flush=True)
print('RESULT Syl16 L=%d: monomials of degree <= %d span all functions on A (mod p): %s' % (L, jdeg, allfull))
