"""Exact check of the averaging constant in Thm 1: for an M-set split into K and M-K, the probability that a uniform
5-subset lies on one side, p(M,K)=(C(K,5)+C(M-K,5))/C(M,5). Claim: min_K p(M,K) >= 1/16 - 2/M for all M>=5."""
from fractions import Fraction as F
from math import comb
worst = None
for M in range(5, 3001):
    p = min(F(comb(K,5)+comb(M-K,5), comb(M,5)) for K in range(0, M+1))
    gap = p - (F(1,16) - F(2, M))
    if worst is None or gap < worst[0]: worst = (gap, M, p)
    assert gap >= 0, (M, p)
print("min_K p(M,K) >= 1/16 - 2/M holds for 5<=M<=3000; tightest margin", float(worst[0]), "at M =", worst[1])
for M in (10, 20, 50, 100, 1000):
    p = min(F(comb(K,5)+comb(M-K,5), comb(M,5)) for K in range(0, M+1))
    print(M, float(p), "1/16 =", 1/16)
