"""Prop 3.1 (invariance/regularity barrier): Q = prod_{i<=m} k(B_i)/sqrt(2b/3) on T={-1,0,1}^d,
disjoint blocks of size b = floor(d/m). Exact: Pr(k(B)=0) = T_b/3^b (central trinomial coefficient, exact integers).
Pr(Q != 0) = (1 - p_b)^m exactly (independent blocks). Also check p_b >= 1/(2 sqrt b) for b <= 5000 (exact)."""
import math
from fractions import Fraction
from common import *

def central_trinomial(b):
    # coefficient of x^b in (1+x+x^2)^b = sum_k b!/(k!k!(b-2k)!)
    s = 0
    c = 1  # k = 0 term: b!/(b!) = 1
    for k in range(0, b // 2 + 1):
        if k == 0:
            c = 1
        else:
            # term_k / term_{k-1} = (b-2k+2)(b-2k+1)/(k^2)
            c = c * (b - 2 * k + 2) * (b - 2 * k + 1) // (k * k)
        s += c
    return s

# lemma check: p_b >= 1/(2 sqrt b), exact rational comparison p_b^2 * 4 b >= 1
bad = []
for b in range(1, 3001):
    T = central_trinomial(b)
    if 25 * b * T * T < 4 * 3 ** (2 * b):  # p_b >= 0.4/sqrt(b)
        bad.append(b)
print("p_b >= 0.4/sqrt(b) fails for b in [1,3000]:", bad[:10], "count", len(bad))

print("n r m' b=floor(d/m') p_b log10 Pr(Q!=0)  [exact p_b; log in 40 digits]  maxInf*d/m'")
for n, r in [(2900, 2900), (4096, 4095), (10**4, 10**4), (65536, 65535), (10**5, 10**5), (10**6, 10**6)]:
    d, w, j, m = wjm(n, r)
    b = d // m
    T = central_trinomial(b)
    p = mpf(T) / mpf(3) ** b
    lg = m * log(1 - p) / log(10)
    print(n, r, m, b, mp.nstr(p, 8), mp.nstr(lg, 8), mp.nstr(mpf(d) / (m * b), 5))
