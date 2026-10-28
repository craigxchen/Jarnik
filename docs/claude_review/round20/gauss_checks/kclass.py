"""Exact: K.gamma_hat and (K + D_PQ).gamma_hat on Mbar_{0,k+2}; and the X_k^# arithmetic
(k-1)(1/2+delta) - 2 (Vojta with D=E, S={inf}) and (k-2)(1/2+delta)-2 (D=0).
K_{Mbar_{0,n}} = sum_S (s(n-s)/(n-1) - 2) D_S (Keel-McKernan/Pandharipande), one term per unordered split.
gamma_hat (threshold/L8.md 6.1): D_PQ -> 1/4 ; D_{P u A}, A nonempty proper subset of [k] -> 2^-k ; else 0."""
from fractions import Fraction as F
from itertools import combinations
def table(k):
    n = k + 2; labels = ['P','Q'] + list(range(1,k+1)); circ = set(range(1,k+1))
    K = F(0); seen = set()
    for r in range(2, n-1):
        for S in combinations(labels, r):
            S = frozenset(S); Sc = frozenset(labels) - S
            key = min(tuple(sorted(map(str,S))), tuple(sorted(map(str,Sc))))
            if key in seen: continue
            seen.add(key)
            s = len(S); c = F(s*(n-s), n-1) - 2
            if S == {'P','Q'} or Sc == {'P','Q'}: g = F(1,4)
            elif ('P' in S) != ('Q' in S): g = F(1, 2**k)
            else: g = F(0)
            K += c*g
    return K, K + F(1,4)
print(" k |  K.gh     | (K+D_PQ).gh | count exp per w = -(K+D).gh*2^(k-1) | X_k: (k-1)/2-2 (D=E) | (k-2)/2-2 (D=0)")
for k in range(3, 10):
    K, KD = table(k)
    print("%2d | %9s | %11s | %6s | %6s | %6s" % (k, K, KD, -KD*2**(k-1), F(k-1,2)-2, F(k-2,2)-2))
# cross-check with direct block count exponent 2^{m-2}(4-m)-1, m=k-1
print("direct count exponent 2^{m-2}(4-m)-1:", [ (k, F(2**(k-1),4)*(4-(k-1))-1) for k in range(3,10)])
