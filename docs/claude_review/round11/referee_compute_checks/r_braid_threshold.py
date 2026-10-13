# ZE-braid vs Prop 5.2 (determinant.md): for the Sylvester-Hadamard profile of order n = 2^k (r = n-1 sign columns)
# and the Theorem 2.1 slab-box (T = n, w = n sqrt(8(n-1) ln 4r)), the k "coordinate" columns S give the subtorus
# delta = sum_{p in S} t_p a_p, contained in Bs(A,m) whenever sum_{p in S} h_A(a_p) < m (Prop 5.2 with kappa = 0).
# On the slice <k,a_p> = <k,beta_p>/2 in [-w/2, w/2] and is an integer, so h_A(a_p) <= 2 floor(w/2).
# Exact m: largest m with 8 binom(m+n-2, n-1) <= (2n+1)^(n-1)  (exact big integers).
import math, sys
from fractions import Fraction
def exact_m(n):
    T = n
    rhs = (2 * T + 1) ** (n - 1)
    lo, hi = int(0.7 * n * (n - 1)), int(0.8 * n * (n - 1))
    assert 8 * math.comb(lo + n - 2, n - 1) <= rhs and 8 * math.comb(hi + n - 2, n - 1) > rhs
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if 8 * math.comb(mid + n - 2, n - 1) <= rhs: lo = mid
        else: hi = mid
    return lo
for k in range(6, 18):
    n = 2 ** k; r = n - 1
    w = n * math.sqrt(8 * (n - 1) * math.log(4 * r))
    hw = 2 * math.floor(w / 2)
    m = exact_m(n) if k <= 16 else None
    m07 = 0.7 * n * (n - 1)
    regime = (n - 1 >= 300 * math.log(4 * r))
    cond07 = k * hw < m07
    line = f"k={k:2d} n={n:6d}: regime(n-1>=300ln4r)={regime}  m>=4w:{(m or m07) >= 4*w}  k*h_A(a_p)<= {k*hw:.4g}  0.7n(n-1)={m07:.4g}"
    if m is not None:
        line += f"  exact m={m}  k*h<m: {k*hw < m}"
    line += f"  (sufficient, via m>=0.7n(n-1)): {cond07}"
    print(line, flush=True)
