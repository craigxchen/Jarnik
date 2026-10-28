"""Exact intersection theory on Mbar_{0,n} (small n) via boundary restriction.
Space = frozenset of labels.  A class (divisor) = dict key -> Fraction, keys ('D', split) or ('P', label)
(split = canonical frozenset side), 'P' = psi.  integrate(L, [c_1..c_d]) = int_{Mbar_{0,L}} c_1...c_d, d=|L|-3.
Restriction to D_T = Mbar_{0,A+h} x Mbar_{0,B+h'}: D_T|D_T = -psi_h - psi_h'; compatible D_S -> boundary of a factor;
incompatible -> 0; psi_l -> psi_l of the factor containing l.  psi_l = sum_{l in T; j,k notin T} D_T.
"""
from fractions import Fraction as Fr
from itertools import combinations
from functools import lru_cache
import itertools
_ctr = itertools.count()
def key(x): return (0, x) if isinstance(x, int) else (1, str(x))
def canon(L, S):
    S = frozenset(S); C = L - S
    m = min(L, key=key)
    return C if m in S else S          # side NOT containing the min label
def all_splits(L):
    L = frozenset(L); m = min(L, key=key); rest = sorted(L - {m}, key=key); out = []
    for r in range(2, len(L) - 1):
        for S in combinations(rest, r):
            S = frozenset(S)
            if len(L) - len(S) >= 2: out.append(S)
    return out
def psi_to_D(L, l):
    others = sorted([x for x in L if x != l], key=key); j, k = others[0], others[1]
    res = {}
    for S in all_splits(L):
        for SS in (S, L - S):
            if l in SS and j not in SS and k not in SS:
                res[('D', canon(L, SS))] = res.get(('D', canon(L, SS)), 0) + 1
    return res
def add(d, k, v):
    if v: d[k] = d.get(k, 0) + v
def restrict(L, T, c, h1, h2):
    A = frozenset(T); B = L - A; X1 = A | {h1}; X2 = B | {h2}
    c1, c2 = {}, {}
    for (typ, x), v in c.items():
        if typ == 'P':
            if x in A: add(c1, ('P', x), v)
            else: add(c2, ('P', x), v)
        else:
            S = x; Sc = L - S
            if S == A or S == B:
                add(c1, ('P', h1), -v); add(c2, ('P', h2), -v); continue
            done = False
            for SS in (S, Sc):
                if SS < A: add(c1, ('D', canon(X1, SS)), v); done = True; break
                if SS < B: add(c2, ('D', canon(X2, SS)), v); done = True; break
            # else incompatible -> 0
    return X1, c1, X2, c2
def freeze(c): return frozenset((k, v) for k, v in c.items() if v)
@lru_cache(maxsize=None)
def _int(L, cs):
    d = len(L) - 3
    assert len(cs) == d
    if d == 0: return Fr(1)
    cs = list(cs); first = dict(cs[0]); rest = [dict(x) for x in cs[1:]]
    # expand psi in first into boundary
    exp = {}
    for (typ, x), v in first.items():
        if typ == 'P':
            for k2, w in psi_to_D(L, x).items(): add(exp, k2, v*w)
        else: add(exp, (typ, x), v)
    tot = Fr(0)
    for (typ, S), v in exp.items():
        if not v: continue
        h1 = ('h', next(_ctr)); h2 = ('h', next(_ctr))
        X1, X2 = None, None
        r1, r2 = [], []
        for c in rest:
            X1, a1, X2, a2 = restrict(L, S, c, h1, h2)
            r1.append(a1); r2.append(a2)
        if X1 is None:
            A = frozenset(S); X1 = A | {h1}; X2 = (L - A) | {h2}
        d1 = len(X1) - 3; d2 = len(X2) - 3
        sub = Fr(0)
        for I in combinations(range(len(rest)), d1):
            J = [j for j in range(len(rest)) if j not in I]
            p1 = integrate(X1, [r1[i] for i in I]); 
            if p1 == 0: continue
            p2 = integrate(X2, [r2[j] for j in J])
            sub += p1 * p2
        tot += v * sub
    return tot
def relabel(L, cs):
    # canonical relabel to ints 0..n-1 to improve caching
    order = sorted(L, key=key); mp = {x: i for i, x in enumerate(order)}
    L2 = frozenset(range(len(order)))
    out = []
    for c in cs:
        c2 = {}
        for (typ, x), v in c.items():
            if typ == 'P': add(c2, ('P', mp[x]), v)
            else: add(c2, ('D', canon(L2, frozenset(mp[y] for y in x))), v)
        out.append(c2)
    return L2, out
def integrate(L, cs):
    L2, cs2 = relabel(frozenset(L), cs)
    if len(L2) - 3 == 0: return Fr(1)
    # order-independent: sort classes for cache
    fz = tuple(sorted((freeze(c) for c in cs2), key=lambda z: sorted(map(str, z))))
    return _int(L2, fz)
def D(L, S): L = frozenset(L); return {('D', canon(L, S)): Fr(1)}
def P(l): return {('P', l): Fr(1)}
