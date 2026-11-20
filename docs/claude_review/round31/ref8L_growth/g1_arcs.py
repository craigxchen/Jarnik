#!/usr/bin/env python3
"""Referee check of VIII.5 incidence bound and VIII.6 affine independence / ratio claims
on genuine circles.  Exact Gaussian integer arithmetic; floats only for sorting by angle.
"""
import math, itertools, random, sys
from fractions import Fraction
from sympy import primerange, factorint
from sympy.matrices import Matrix

def gprime(p):
    # x^2+y^2=p, x>y>0
    for y in range(1, int(math.isqrt(p)) + 1):
        x2 = p - y*y
        x = math.isqrt(x2)
        if x*x == x2:
            return (max(x, y), min(x, y))
    raise ValueError

def gmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])

def gpow(a, e):
    r = (1, 0)
    for _ in range(e):
        r = gmul(r, a)
    return r

def gconj(a):
    return (a[0], -a[1])

def gnorm(a):
    return a[0]*a[0]+a[1]*a[1]

def gdivmod_exact(a, b):
    # a/b if exact else None
    n = gnorm(b)
    c = gmul(a, gconj(b))
    if c[0] % n or c[1] % n:
        return None
    return (c[0]//n, c[1]//n)

def ggcd(a, b):
    while b != (0, 0):
        n = gnorm(b)
        c = gmul(a, gconj(b))
        q = (round(Fraction(c[0], n)), round(Fraction(c[1], n)))
        qb = gmul(q, b)
        a, b = b, (a[0]-qb[0], a[1]-qb[1])
    return a

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def points(fac):
    """fac: dict p->e of split primes. returns list of (z, valuation vector a) for all z with |z|^2=N,
    primitive in the sense gcd(z, conj z) may be nontrivial (we keep all)."""
    ps = sorted(fac)
    pis = {p: gprime(p) for p in ps}
    res = []
    for avec in itertools.product(*[range(fac[p]+1) for p in ps]):
        z = (1, 0)
        for p, a in zip(ps, avec):
            z = gmul(z, gmul(gpow(pis[p], a), gpow(gconj(pis[p]), fac[p]-a)))
        for u in UNITS:
            res.append((gmul(u, z), avec))
    return ps, res

def arg(z):
    return math.atan2(z[1], z[0]) % (2*math.pi)

def clusters(pts, width):
    pts = sorted(pts, key=lambda t: arg(t[0]))
    angs = [arg(t[0]) for t in pts]
    n = len(pts)
    out = []
    j = 0
    for i in range(n):
        cl = [pts[i]]
        k = 1
        while k < n and ((angs[(i+k) % n] - angs[i]) % (2*math.pi)) <= width:
            cl.append(pts[(i+k) % n]); k += 1
        out.append(cl)
    return out

def reduce_ratio(z, w):
    g = ggcd(z, w)
    return gdivmod_exact(z, g), gdivmod_exact(w, g)

def normalize_unit(A, B):
    # canonical representative of A/B up to common unit
    best = None
    for u in UNITS:
        c = (gmul(u, A), gmul(u, B))
        if best is None or c > best:
            best = c
    return best

def main():
    random.seed(1)
    split = [p for p in primerange(5, 400) if p % 4 == 1]
    stats = dict(circles=0, maxpairs_ratio=0, ratio_violations=0, hyper_checked=0, hyper_viol=0,
                 ai_clusters=0, ai_fail=0, Dviol=0, ratio_dup=0, ratio_checked=0, height_viol=0, maxm=0)
    worst_ratio = {}
    for trial in range(int(sys.argv[1]) if len(sys.argv) > 1 else 200):
        k = random.randint(3, 7)
        ps = random.sample(split[:25], k)
        fac = {p: random.choice([1, 1, 1, 2, 3]) for p in ps}
        N = 1
        for p, e in fac.items():
            N *= p**e
        ps, pts = points(fac)
        stats['circles'] += 1
        W = math.log(N)
        for C in (1.0, 2.0, 3.0):
            width = C * N**-0.25
            # ratio bound: all pairs in some arc of width `width`
            for cl in clusters(pts, width):
                ratio_count = {}
                seen_pairs = set()
                for a, b in itertools.permutations(range(len(cl)), 2):
                    z, w = cl[a][0], cl[b][0]
                    if (z, w) in seen_pairs:
                        continue
                    seen_pairs.add((z, w))
                    A, B = reduce_ratio(z, w)
                    key = normalize_unit(A, B)
                    ratio_count[key] = ratio_count.get(key, 0) + 1
                if ratio_count:
                    mx = max(ratio_count.values())
                    stats['maxpairs_ratio'] = max(stats['maxpairs_ratio'], mx)
                    worst_ratio[C] = max(worst_ratio.get(C, 0), mx)
                    stats['maxM_C%g' % C] = max(stats.get('maxM_C%g' % C, 0), len(cl))
                    if mx > 1 + C*C:
                        stats['ratio_violations'] += 1
        # VIII.6: clusters of width < N^{-1/4}
        width = 0.999 * N**-0.25
        for cl in clusters(pts, width):
            m = len(cl)
            if m < 2:
                continue
            stats['maxm'] = max(stats['maxm'], m)
            stats['ai_clusters'] += 1
            rows = [list(t[1]) + [1] for t in cl]
            if Matrix(rows).rank() < m:
                stats['ai_fail'] += 1
            for (z, a), (w, b) in itertools.combinations(cl, 2):
                D = sum(abs(x-y)*math.log(p) for x, y, p in zip(a, b, ps))
                if not D > W/2:
                    stats['Dviol'] += 1
            # ratios for l=1 (and l=2 if m>=4)
            for l in (1, 2):
                if m < 2*l:
                    continue
                seen = {}
                for I in itertools.combinations(range(m), l):
                    rest = [x for x in range(m) if x not in I]
                    for J in itertools.combinations(rest, l):
                        num = (1, 0); den = (1, 0)
                        for i in I: num = gmul(num, cl[i][0])
                        for j in J: den = gmul(den, cl[j][0])
                        A, B = reduce_ratio(num, den)
                        key = normalize_unit(A, B)
                        stats['ratio_checked'] += 1
                        if key in seen:
                            stats['ratio_dup'] += 1
                        seen[key] = 1
                        H = gnorm(B)
                        u = complex(*A) / complex(*B)
                        if not (H <= N**l and H >= math.sqrt(N)/l**2 * 0.999999 and abs(u-1) <= l*N**-0.25*1.000001):
                            stats['height_viol'] += 1
    print(stats)
    print('max pairs per ratio by C:', worst_ratio)

if __name__ == '__main__':
    main()
