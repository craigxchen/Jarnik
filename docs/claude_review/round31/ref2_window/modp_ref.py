"""Referee's independent finite-field test of the degree-one K5 moment system.

Field F_{p^2} = F_p[i], p = 3 mod 4, conjugation = Frobenius (a+bi -> a-bi).
Unknowns: ten roots r_e (e an edge of K5), all twenty r_e, conj(r_e) distinct (so nonreal).
Equations: for k = 1..KMAX the quantity S_v(k) = sum_{e ni v} Im(r_e^k) is the same for all
five vertices v.  (Degree-one case of (II.8), KMAX = 6d-1 = 5.)

Algorithm (different from Astra's): build all "stars" = 4 distinct conjugate pairs with an
orientation; bucket by moment vector; inside each bucket find all 5-sets of stars that form an
edge-labelled K5: pairwise intersections of size exactly one, the 10 intersection elements
distinct, and every element of the union in exactly 2 stars (forced by counting), plus the
conjugacy condition (the 10 used conjugate pairs distinct), unless --loose, which only asks that
the 10 oriented roots be distinct and each star avoid its own conjugates.
Each solution is then re-verified directly from the edge labelling.
"""
import sys
from itertools import combinations, product
from collections import defaultdict

p = int(sys.argv[1]); KMAX = int(sys.argv[2]); loose = '--loose' in sys.argv
assert p % 4 == 3

def mul(z, w):
    return ((z[0]*w[0]-z[1]*w[1]) % p, (z[0]*w[1]+z[1]*w[0]) % p)

elems = [(a, b) for a in range(p) for b in range(1, p)]          # nonreal elements
pairs = sorted({tuple(sorted([(a, b), (a, (-b) % p)])) for (a, b) in elems})
assert len(pairs) == p*(p-1)//2
def immom(z):
    out = []; w = (1, 0)
    for k in range(KMAX):
        w = mul(w, z); out.append(w[1])
    return out
mom = {z: immom(z) for z in elems}
conj = {z: (z[0], (-z[1]) % p) for z in elems}

buckets = defaultdict(list)
nst = 0
for quad in combinations(range(len(pairs)), 4):
    for sg in product((0, 1), repeat=4):
        star = frozenset(pairs[q][s] for q, s in zip(quad, sg))
        key = tuple(sum(mom[z][k] for z in star) % p for k in range(KMAX))
        buckets[key].append(star); nst += 1
print('p', p, 'KMAX', KMAX, 'stars', nst, 'classes', len(buckets),
      'max class', max(len(b) for b in buckets.values()), flush=True)

def verify(stars):
    # stars: list of 5 frozensets; derive labelling and check everything from scratch
    lab = {}
    for a, b in combinations(range(5), 2):
        I = stars[a] & stars[b]
        assert len(I) == 1
        lab[(a, b)] = next(iter(I))
    roots = list(lab.values())
    assert len(set(roots)) == 10
    if not loose:
        allz = roots + [conj[z] for z in roots]
        assert len(set(allz)) == 20
    for v in range(5):
        assert set(lab[e] for e in lab if v in e) == set(stars[v])
    S = [tuple(sum(mom[lab[e]][k] for e in lab if v in e) % p for k in range(KMAX)) for v in range(5)]
    assert len(set(S)) == 1
    return True

found = 0
for key, B in buckets.items():
    n = len(B)
    if n < 5: continue
    ok = [[False]*n for _ in range(n)]
    for a in range(n):
        for b in range(a+1, n):
            if len(B[a] & B[b]) != 1: continue
            if not loose:
                # no element of one star may be the conjugate of an element of the other
                if any(conj[z] in B[b] for z in B[a]): continue
            ok[a][b] = ok[b][a] = True
    def rec(chosen, start):
        global found
        if len(chosen) == 5:
            st = [B[c] for c in chosen]
            inter = [next(iter(st[a] & st[b])) for a, b in combinations(range(5), 2)]
            if len(set(inter)) != 10: return
            if not loose:
                pr = [frozenset([z, conj[z]]) for z in inter]
                if len(set(pr)) != 10: return
            verify(st); found += 1
            if found <= 3: print('solution', [sorted(s) for s in st], flush=True)
            return
        for c in range(start, n):
            if all(ok[c][d] for d in chosen):
                rec(chosen+[c], c+1)
    rec([], 0)
print('unordered solutions (sets of 5 stars):', found, 'loose' if loose else 'strict')
