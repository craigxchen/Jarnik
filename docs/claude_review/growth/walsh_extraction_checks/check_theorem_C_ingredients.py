"""Exact checks of the ingredients of Theorems B'' and C (walsh-extraction.md).
1. Density lemma (item 326): constructive doubling finds a K-point affine flat inside any A of density
   >= 2 (K/M)^{1/K}; checked on random sets.
2. Row-deleted pigeonhole certificate: if (Hv)|_F = (Hv')|_F then lambda = H(v-v')/M is supported on U,
   has zero sum, lies in (1/M) Z^M and has column images exactly v-v' (all exact rationals).
3. Saturation: in a fully flipped profile (one flip per row, distinct columns, an unflipped copy of every
   label) every rational lambda with integral column images has 2*lambda integral; in a partially flipped
   Walsh profile an affine sub-cube Y of unflipped rows carries the fractional characters chi_a|_Y / |Y|.
"""
import random, sys, os
from fractions import Fraction as Fr
random.seed(11)

def walsh(t):
    M = 2 ** t
    return [[(-1) ** bin(x & a).count("1") for a in range(M)] for x in range(M)]

def paley(q):
    sq = set((a * a) % q for a in range(1, q))
    ch = lambda a: 0 if a % q == 0 else (1 if (a % q) in sq else -1)
    M = q + 1
    H = [[1] * M for _ in range(M)]
    for x in range(1, M):
        for j in range(1, M):
            H[x][j] = -1 if x == j else -ch(x - j)
    return H

# 1. density lemma
def find_flat(A, t, K):
    M = 2 ** t
    cosets = [frozenset([a]) for a in A]       # cosets of V={0}
    V = [0]
    h = 1
    while h < K:
        n = len(cosets)
        if n < 2: return None
        reps = [min(c) for c in cosets]
        # quotient by V: canonical rep = min of coset; direction between cosets P,Q: (rep P xor rep Q) reduced mod V
        def red(d):
            return min(d ^ v for v in V)
        cnt = {}
        for i in range(n):
            for j in range(n):
                if i != j:
                    d = red(reps[i] ^ reps[j]); cnt[d] = cnt.get(d, 0) + 1
        e = max(cnt, key=cnt.get)
        setc = set(cosets); used = set(); new = []
        for c in cosets:
            if c in used: continue
            partner = frozenset(x ^ e for x in c)
            if partner in setc and partner not in used and partner != c:
                used.add(c); used.add(partner); new.append(c | partner)
        V = sorted(set(V) | set(v ^ e for v in V))
        cosets = new; h *= 2
    return cosets[0] if cosets else None

ok1 = 0
for t in range(6, 12):
    M = 2 ** t
    for K in (2, 4, 8):
        alpha = 2 * (K / M) ** (1.0 / K)
        if alpha > 0.95: continue
        for trial in range(5):
            A = [x for x in range(M) if random.random() < min(1.0, alpha + 0.02)]
            if len(A) < alpha * M: continue
            F = find_flat(set(A), t, K)
            assert F is not None and len(F) == K and F <= set(A)
            base = min(F); D = set(x ^ base for x in F)
            assert all((u ^ v) in D for u in D for v in D)     # linear subspace => affine flat
            ok1 += 1
print("1. density lemma: affine flats found in", ok1, "random dense sets (all at the stated threshold)")

# 2. row-deleted certificate
ok2 = 0
for name, H in [("walsh32", walsh(5)), ("paley44", paley(43)), ("walsh64", walsh(6))]:
    M = len(H)
    for trial in range(30):
        f = random.randint(1, 4)
        Fr_rows = random.sample(range(M), f)
        U = [x for x in range(M) if x not in Fr_rows]
        L = random.randint(2, 6)
        seen = {}
        found = None
        for it in range(20000):
            v = [0] * M
            for a in random.sample(range(1, M), L): v[a] = 1
            key = tuple(sum(H[x][a] * v[a] for a in range(1, M)) for x in Fr_rows)
            if key in seen and seen[key] != v:
                found = (seen[key], v); break
            seen[key] = v
        if not found: continue
        v, v2 = found
        n = [v[a] - v2[a] for a in range(M)]
        lam = [Fr(sum(H[x][a] * n[a] for a in range(1, M)), M) for x in range(M)]
        assert all(lam[x] == 0 for x in Fr_rows)
        assert sum(lam) == 0
        assert all((lam[x] * M).denominator == 1 for x in range(M))
        for a in range(1, M):
            assert sum(lam[x] * H[x][a] for x in U) == n[a]
        ok2 += 1
print("2. row-deleted pigeonhole certificates verified exactly:", ok2)

# 3. saturation
def flipped(H, b, flip_rows):
    M = len(H); cols = []; lab = []
    for a in range(1, M):
        for c in range(b): cols.append([H[x][a] for x in range(M)]); lab.append(a)
    order = list(flip_rows); random.shuffle(order)
    for idx, x in enumerate(order):
        a = 1 + idx % (M - 1); j = (a - 1) * b + idx // (M - 1)
        cols[j][x] *= -1
    return cols, lab
ok3 = 0
for name, H in [("walsh16", walsh(4)), ("paley12", paley(11)), ("paley20", paley(19))]:
    M = len(H)
    cols, lab = flipped(H, 5, range(M))
    # every row: a flipped column and an unflipped copy of its label differ exactly by -+2 e_x
    for x in range(M):
        js = [j for j in range(len(cols)) if cols[j][x] != H[x][lab[j]]]
        assert len(js) == 1
        j = js[0]
        sib = [k for k in range(len(cols)) if lab[k] == lab[j] and k != j and all(cols[k][y] == H[y][lab[k]] for y in range(M))]
        assert sib
        d = [cols[j][y] - cols[sib[0]][y] for y in range(M)]
        assert all(d[y] == 0 for y in range(M) if y != x) and abs(d[x]) == 2
        ok3 += 1
print("3a. fully flipped: every row has a flipped/unflipped sibling pair differing by 2e_x (so 2*lambda_x integral):", ok3, "rows")
# partially flipped Walsh: fractional characters on an unflipped affine sub-cube
H = walsh(6); M = 64
Y = [x for x in range(M) if (x & 0b110000) == 0b010000]          # affine 4-cube (16 rows)
flip_rows = [x for x in range(M) if x not in Y]
cols, lab = flipped(H, 5, flip_rows)
ok4 = 0
for a in range(1, M):
    # restricted label of a on Y: Y = x0 + span(e0..e3); chi_a on Y
    lam = [Fr(0)] * M
    for x in Y: lam[x] = Fr(H[x][a], len(Y))
    if sum(lam) != 0: continue          # label constant on Y
    imgs = [sum(lam[x] * cols[j][x] for x in range(M)) for j in range(len(cols))]
    assert all(im.denominator == 1 for im in imgs)
    assert any(im != 0 for im in imgs)
    ok4 += 1
print("3b. partially flipped Walsh-64 with an unflipped 16-row affine cube: %d fractional (1/16) characters with integral images" % ok4)
print("ALL PASS")
