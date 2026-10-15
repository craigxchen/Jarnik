# Referee check: with GENUINE labels at all odd q <= X (X >= a_j+b_j), are columns pinned at norm 1
# using all four units s=0..3 (not only S(theta))?  Solution set of (1.2) for e_j is [0,pi) minus four
# open arcs centred at s*pi/4 with radii arcsin(nu_s m_s/sqrt p), m_s = odd, N-free part of D_s.
# Column pinned iff every gap between consecutive arcs has length 0.  Gaps are computed in mpmath at
# 60 digits; "zero" means |gap| < 1e-40 (an exact zero is the Gaussian identity of a tight pair).
from sympy import primerange
from math import isqrt
import mpmath as mp
mp.mp.dps = 60

def gauss(p):
    for b in range(2, isqrt(p)+1, 2):
        a2 = p - b*b; a = isqrt(a2)
        if a*a == a2 and a % 2: return a, b

def nfree(n, ps):
    n = abs(n)
    while n % 2 == 0: n //= 2
    for pl in ps:
        while n % pl == 0: n //= pl
    return n

def gaps(p, ps):
    # complement of the union of open arcs (c_s - r_s, c_s + r_s) mod pi; returns list of complement lengths
    a, b = gauss(p)
    D = [b, a-b, a, a+b]
    nu = [1, 1/mp.sqrt(2), 1, 1/mp.sqrt(2)]
    ivs = []
    for s in range(4):
        r = mp.asin(nu[s]*nfree(D[s], ps)/mp.sqrt(p)); c = s*mp.pi/4
        lo, hi = c - r, c + r
        for sh in (-mp.pi, 0, mp.pi):
            L, H = max(lo+sh, mp.mpf(0)), min(hi+sh, mp.pi)
            if L < H: ivs.append((L, H))
    ivs.sort()
    comp = []; cur = mp.mpf(0)
    for L, H in ivs:
        if L > cur: comp.append(L - cur)
        cur = max(cur, H)
    if cur < mp.pi: comp.append(mp.pi - cur)
    # glue wrap-around component
    return comp if comp else [mp.mpf(0)]

def unpinned(ps):
    out = []
    for p in ps:
        g = gaps(p, ps)
        if max(g) > mp.mpf(10)**-40:
            out.append((p, float(max(g))))
    return out

split = [p for p in primerange(3, 10**6) if p % 4 == 1]
print("first r split primes, genuine labels at every odd q <= X, X >= max(a_j+b_j):")
for r in [10, 20, 50, 100, 200, 500]:
    ps = split[:r]
    u = unpinned(ps)
    print(f"  r={r}: unpinned columns at norm 1 (all four units) = {len(u)} / {r}; examples {[(p, round(g,4)) for p,g in u[:4]]}")
print("sanity: primes in [10^4,10^5] (random 200): ", end="")
import random; random.seed(3)
ps = random.sample([p for p in split if 10**4 < p < 10**5], 200)
print(len(unpinned(ps)), "/ 200 unpinned")

print("round-10 L4K prime set {5,13,17,29} (flexible-top instance):")
ps = [5, 13, 17, 29]
for p in ps:
    a, b = gauss(p)
    print(f"  p={p}: a={a}, b={b}, D=(b,a-b,a,a+b)=({b},{a-b},{a},{a+b}); free-set lengths mod pi = {[round(float(g),4) for g in gaps(p, ps)]}")
print("  unpinned:", unpinned(ps))
