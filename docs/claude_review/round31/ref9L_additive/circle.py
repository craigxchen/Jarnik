#!/usr/bin/env python3
"""Referee check of IX.3 (circle) on genuine circles x^2+y^2=N, N composed of
split primes. Points z = eps * prod pi_p^{a_p} conj(pi_p)^{e_p-a_p} (exact
Gaussian integers, exponent rows and unit recorded).
For each window of angular width theta = C*N^{-1/4} (arc length C*N^{1/4}=C sqrt R):
  * M = number of points, normalised exponents a-lo, e'_p = hi-lo, N' = prod p^{e'_p};
  * unary signatures in F_2^{Omega(N')}, D = |A+A|; record M and sqrt(D);
  * core lemma test: inside one unit class and a subarc of angular diameter
    c*N'^{-1/4}, unordered pairs have distinct signature differences
    (asserted for c<2 by my recomputed gap 2e^{-w/2}; Astra asserts c<1/10).
    We test c = 1.99 and also record violations at larger c (2.5, 4, 8).
  * gap check |alpha(v) -+ 1| >= 2 e^{-w(v)/2} on all pairs in a unit class.
"""
import sys, math, itertools
from sympy import factorint, primerange
from collections import Counter

def gauss_prime(p):
    for a in range(1, math.isqrt(p) + 1):
        b2 = p - a * a
        b = math.isqrt(b2)
        if b * b == b2:
            return (a, b)
    raise ValueError

def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

def gpow(x, k):
    r = (1, 0)
    for _ in range(k):
        r = gmul(r, x)
    return r

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def points(fac):
    ps = sorted(fac)
    pis = {p: gauss_prime(p) for p in ps}
    pw = {p: [gmul(gpow(pis[p], a), gpow((pis[p][0], -pis[p][1]), fac[p] - a)) for a in range(fac[p] + 1)] for p in ps}
    pts = []
    for row in itertools.product(*[range(fac[p] + 1) for p in ps]):
        z = (1, 0)
        for p, a in zip(ps, row):
            z = gmul(z, pw[p][a])
        for u, eps in enumerate(UNITS):
            w = gmul(eps, z)
            pts.append((math.atan2(w[1], w[0]) % (2 * math.pi), u, row, w))
    pts.sort()
    return ps, pts

def analyse(fac, C, cs=(1.99, 2.5, 4.0, 8.0)):
    ps, pts = points(fac)
    N = 1
    for p in ps: N *= p ** fac[p]
    theta = C * N ** -0.25
    n = len(pts)
    out0 = dict(maxM=0, maxratio=0.0, viol={c: 0 for c in cs}, gapviol=0, windows=0, recs=[])
    if theta >= 1.5:
        return N, out0
    ext = pts + [(t + 2 * math.pi, u, r, w) for (t, u, r, w) in pts]
    out = dict(maxM=0, maxratio=0.0, viol={c: 0 for c in cs}, gapviol=0, windows=0, recs=[])
    j = 0
    for i in range(n):
        while j < i + n and ext[j][0] - ext[i][0] <= theta:
            j += 1
        W = ext[i:j]
        M = len(W)
        if M < 2: continue
        out['windows'] += 1
        lo = [min(r[k] for _, _, r, _ in W) for k in range(len(ps))]
        hi = [max(r[k] for _, _, r, _ in W) for k in range(len(ps))]
        Np = 1
        for k, p in enumerate(ps): Np *= p ** (hi[k] - lo[k])
        # unary signature
        def sig(r):
            v = 0; off = 0
            for k in range(len(ps)):
                e = hi[k] - lo[k]
                a = r[k] - lo[k]
                v |= ((1 << a) - 1) << off
                off += e
            return v
        A = [sig(r) for _, _, r, _ in W]
        assert len(set(A)) == M
        D = len(set(x ^ y for x in A for y in A))
        out['maxM'] = max(out['maxM'], M)
        rr = M / math.sqrt(D)
        if rr > out['maxratio']: out['maxratio'] = rr
        out['recs'].append((M, D))
        # core lemma and gap, within unit classes
        for u in range(4):
            Wu = [(t, r, a) for (t, uu, r, _), a in zip(W, A) if uu == u]
            for c in cs:
                ang = c * Np ** -0.25
                # all subarcs starting at each point of Wu
                for s0 in range(len(Wu)):
                    sub = [x for x in Wu[s0:] if x[0] - Wu[s0][0] <= ang]
                    if len(sub) < 3: continue
                    d = Counter(x[2] ^ y[2] for x, y in itertools.combinations(sub, 2))
                    if max(d.values()) > 1:
                        out['viol'][c] += 1
            for (t1, r1, _), (t2, r2, _) in itertools.combinations(Wu, 2):
                v = [x - y for x, y in zip(r1, r2)]
                w = sum(abs(x) * math.log(p) for x, p in zip(v, ps))
                dang = abs(t1 - t2)
                # |alpha-1| = 2 sin(dang/2) ; must be >= 2 e^{-w/2}
                if 2 * math.sin(dang / 2) < 2 * math.exp(-w / 2) * (1 - 1e-9):
                    out['gapviol'] += 1
    return N, out

if __name__ == '__main__':
    mode = sys.argv[1]
    Cs = [float(x) for x in sys.argv[2].split(',')]
    if mode == 'struct':
        split = [p for p in primerange(5, 400) if p % 4 == 1]
        facs = []
        for k in range(3, 11):
            facs.append({p: 1 for p in split[:k]})
        for k in range(2, 7):
            facs.append({p: 2 for p in split[:k]})
        for k in range(2, 5):
            facs.append({p: 3 for p in split[:k]})
        facs.append({5: 6, 13: 3, 17: 2, 29: 1, 37: 1})
        facs.append({5: 4, 13: 4, 17: 2, 29: 2})
        facs.append({p: 1 for p in split[10:19]})
    else:
        lo, hi = int(sys.argv[3]), int(sys.argv[4])
        facs = []
        for N in range(lo, hi + 1):
            f = factorint(N)
            if all(p % 4 == 1 for p in f) and len(f) >= 3:
                facs.append(f)
    for C in Cs:
        tot = Counter(); best = (0, None); bestM = (0, None); hist = Counter()
        for fac in facs:
            N, o = analyse(fac, C)
            for c, v in o['viol'].items(): tot['viol_c=%g' % c] += v
            tot['gapviol'] += o['gapviol']; tot['windows'] += o['windows']
            if o['maxratio'] > best[0]: best = (o['maxratio'], fac)
            if o['maxM'] > bestM[0]: bestM = (o['maxM'], fac)
            for M, D in o['recs']: hist[(M, D)] += 1
            if mode == 'struct':
                mr = max(o['recs'], key=lambda x: x[0]) if o['recs'] else None
                print('  C=%g N=%s maxM=%d (M,D) at maxM=%s max M/sqrtD=%.3f viol=%s gap=%d' % (
                    C, fac, o['maxM'], mr, o['maxratio'], o['viol'], o['gapviol']), flush=True)
        print('C=%g' % C, dict(tot), 'max M/sqrt(D)=%.3f at %s' % best, 'max M=%d at %s' % bestM)
        print('   (M,D) histogram top:', sorted(hist.items(), key=lambda x: (-x[0][0], x[0][1]))[:12], flush=True)
