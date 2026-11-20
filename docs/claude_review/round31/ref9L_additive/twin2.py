#!/usr/bin/env python3
"""DFS version: enumerate maximal admissible configurations S in [0,L]
(at most two classes mod each p^f || q) with |S| >= thr, by branching on the
class pair at each prime power among residues actually present. Same checks
as twin.py (signature injectivity, affine 3-spaces, difference multiplicities,
sumset bound, B/sqrt(D))."""
import sys, math
from collections import Counter
from sympy import factorint
from twin import has_affine3

def configs(q, L, thr, cap=200000):
    pps = sorted((p ** f for p, f in factorint(q).items()), reverse=True)
    out = []
    def rec(i, cand, ch):
        if len(out) >= cap: return
        if i == len(pps):
            out.append((tuple(cand), tuple(ch))); return
        m = pps[i]
        cnt = Counter(y % m for y in cand)
        rs = sorted(cnt, key=lambda r: -cnt[r])
        if len(rs) <= 1:
            rs2 = [(rs[0], None)] if rs else []
        else:
            rs2 = [(a, b) for ia, a in enumerate(rs) for b in rs[ia+1:]]
        for a, b in rs2:
            c = cnt[a] + (cnt[b] if b is not None else 0)
            if c < thr: continue
            rec(i + 1, [y for y in cand if y % m in (a, b)], ch + [(m, a, b)])
    rec(0, list(range(L + 1)), [])
    return pps, out

def check(q, Cn, Cd, thr, cap):
    L = math.isqrt(Cn * Cn * q // (Cd * Cd)); Cf = Cn / Cd
    assert L < q
    pps, out = configs(q, L, thr, cap)
    mb = (1 + Cf * Cf) * (1 + 2 * Cf * Cf)
    res = Counter(); worst = dict(B=0, mult=0, ratio=0, Bex=None)
    for S, ch in out:
        if len(S) < 2: continue
        A = []
        for y in S:
            v = 0
            for t, (m, a, b) in enumerate(ch):
                if b is not None and y % m == b: v |= 1 << t
            A.append(v)
        assert len(set(A)) == len(A)
        B = len(A); res['conf'] += 1
        d = Counter(A[i] ^ A[j] for i in range(B) for j in range(i+1, B))
        D = len(set(a ^ b for a in A for b in A))
        mm = max(d.values())
        if mm > mb + 1e-9: res['viol_mult'] += 1
        for h, c in d.items():
            H = 1
            for t, (m, a, b) in enumerate(ch):
                if not (h >> t) & 1: H *= m
            K = q // H
            if c > (2 * L // K + 1) * ((L // H) ** 2 // K + 1): res['viol_sharp'] += 1
        if B * (B - 1) > (D - 1) * (1 + math.sqrt(1 + 8 * (D - 1))) + 1e-9: res['viol_sum'] += 1
        if B >= 8 and has_affine3(A):
            res['aff3'] += 1
            if q > (4 * Cf) ** 6: res['aff3_large'] += 1
        if B > worst['B']: worst['B'] = B; worst['Bex'] = (S, ch, D)
        worst['mult'] = max(worst['mult'], mm)
        worst['ratio'] = max(worst['ratio'], B / math.sqrt(D))
    return L, len(pps), res, worst

if __name__ == '__main__':
    Cn, Cd, thr, cap = map(int, sys.argv[1:5])
    for q in map(int, sys.argv[5:]):
        L, s, res, w = check(q, Cn, Cd, thr, cap)
        print('q=%d L=%d s=%d C=%d/%d thr=%d' % (q, L, s, Cn, Cd, thr), dict(res),
              'maxB=%d maxmult=%d max B/sqrtD=%.3f' % (w['B'], w['mult'], w['ratio']),
              'bound=%.1f' % ((1+(Cn/Cd)**2)*(1+2*(Cn/Cd)**2)), flush=True)
        if w['Bex'] and w['B'] >= 6: print('   example S=', w['Bex'][0], 'D=', w['Bex'][2])
