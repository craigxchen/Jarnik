#!/usr/bin/env python3
"""Referee check of IX.1 and IX.3 (twin) on small moduli.
For q and C, L = floor(C*sqrt(q)). By translation invariance every admissible
configuration is a translate of a set S in [0,L] meeting at most two classes
mod every exact prime power p^f || q. Maximal S are given by a class pair at
each prime power. We enumerate all class pairs (or sample when too many).
Checks per configuration:
  (a) signature injective;
  (b) no affine 3-space in A (only asserted for q > (4C)^6) -- we record all;
  (c) multiplicity of each signature difference <= (1+C^2)(1+2C^2) and <= the
      sharper exact bound (1+2L/K)(1+L^2/(H^2 K)) floor-counted;
  (d) B(B-1) <= (D-1)(1+sqrt(1+8(D-1)));
  (e) record max B, D, B/sqrt(D).
"""
import sys, math, random, itertools
from sympy import factorint
from collections import Counter

def isqrt_floor_mul(C, q):
    # floor(C*sqrt(q)) with C rational given as (num,den)
    num, den = C
    return math.isqrt(num * num * q // (den * den))

def has_affine3(A):
    S = set(A)
    L = sorted(A)
    for a in L:
        T = [a ^ b for b in L if b > a]  # space with min element a: all others > a
        Ts = set(T)
        n = len(T)
        for i in range(n):
            u1 = T[i]
            for j in range(i + 1, n):
                u2 = T[j]
                u12 = u1 ^ u2
                if u12 not in Ts:
                    continue
                for k in range(j + 1, n):
                    u3 = T[k]
                    if u3 == u12:
                        continue
                    if (u3 ^ u1) in Ts and (u3 ^ u2) in Ts and (u3 ^ u12) in Ts:
                        return (a, u1, u2, u3)
    return None

def analyse(q, C, maxchoices=200000, rng=random.Random(1)):
    fac = factorint(q)
    pps = sorted(p ** f for p, f in fac.items())
    L = isqrt_floor_mul(C, q)
    ys = list(range(L + 1))
    pairs_per = [list(itertools.combinations(range(m), 2)) for m in pps]
    total = 1
    for pl in pairs_per:
        total *= len(pl)
    if total <= maxchoices:
        choices = itertools.product(*pairs_per)
        exhaustive = True
    else:
        choices = (tuple(rng.choice(pl) for pl in pairs_per) for _ in range(maxchoices))
        exhaustive = False
    Cf = C[0] / C[1]
    multbound = (1 + Cf * Cf) * (1 + 2 * Cf * Cf)
    stats = dict(q=q, L=L, s=len(pps), exhaustive=exhaustive, nconf=0, maxB=0, maxmult=0,
                 viol_mult=0, viol_sharp=0, viol_sum=0, affine3=0, affine3_large=0,
                 maxratio=0.0, ex_aff=None)
    seen = set()
    for ch in choices:
        S = [y for y in ys if all((y % m) in c for m, c in zip(pps, ch))]
        key = tuple(S)
        if len(S) < 2 or key in seen:
            continue
        seen.add(key)
        stats['nconf'] += 1
        sig = {}
        for y in S:
            v = 0
            for t, (m, c) in enumerate(zip(pps, ch)):
                if y % m == c[1]:
                    v |= 1 << t
            sig[y] = v
        A = list(sig.values())
        assert len(set(A)) == len(A), "signature not injective"
        B = len(A)
        diffs = Counter()
        for i in range(B):
            for j in range(i + 1, B):
                diffs[A[i] ^ A[j]] += 1
        D = len(set(a ^ b for a in A for b in A))
        mm = max(diffs.values())
        stats['maxmult'] = max(stats['maxmult'], mm)
        if mm > multbound + 1e-9:
            stats['viol_mult'] += 1
        # sharper exact bound per h
        for h, cnt in diffs.items():
            H = 1
            for t, m in enumerate(pps):
                if not (h >> t) & 1:
                    H *= m
            K = q // H
            b1 = 2 * L // K + 1
            b2 = (L // H) ** 2 // K + 1
            if cnt > b1 * b2:
                stats['viol_sharp'] += 1
        if B * (B - 1) > (D - 1) * (1 + math.sqrt(1 + 8 * (D - 1))) + 1e-9:
            stats['viol_sum'] += 1
        if B >= 8:
            af = has_affine3(A)
            if af:
                stats['affine3'] += 1
                if q > (4 * Cf) ** 6:
                    stats['affine3_large'] += 1
                if stats['ex_aff'] is None:
                    stats['ex_aff'] = (ch, S)
        stats['maxB'] = max(stats['maxB'], B)
        stats['maxratio'] = max(stats['maxratio'], B / math.sqrt(D))
    return stats

if __name__ == '__main__':
    Cnum, Cden = int(sys.argv[1]), int(sys.argv[2])
    qlo, qhi = int(sys.argv[3]), int(sys.argv[4])
    maxch = int(sys.argv[5]) if len(sys.argv) > 5 else 20000
    agg = Counter()
    worst = {}
    for q in range(qlo, qhi + 1):
        st = analyse(q, (Cnum, Cden), maxch)
        for k in ('nconf', 'viol_mult', 'viol_sharp', 'viol_sum', 'affine3', 'affine3_large'):
            agg[k] += st[k]
        for k in ('maxB', 'maxmult', 'maxratio'):
            if st[k] > worst.get(k, (0, None))[0]:
                worst[k] = (st[k], q)
        if st['affine3'] and 'aff_example' not in worst:
            worst['aff_example'] = (q, st['L'], st['ex_aff'])
    print('C=%d/%d q in [%d,%d]' % (Cnum, Cden, qlo, qhi), dict(agg))
    print('worst', worst)
