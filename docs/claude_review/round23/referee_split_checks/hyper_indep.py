"""Check Thm 1 step (ii): p(M,k)=(C(k,5)+C(M-k,5))/C(M,5); min over k vs (1/16)(1-15/(M-4)) (M>=9) and 1/16-2/M.
Also k-subset version: min_k p_k(M) >= 2^{1-k}(1-O(1/M)) for k=6,7 (Remark 2)."""
from fractions import Fraction as Fr
from math import comb
bad1=bad2=0
for M in range(5,3001):
    mn=min(Fr(comb(k,5)+comb(M-k,5),comb(M,5)) for k in range(M+1))
    if M>=9 and mn < Fr(1,16)*(1-Fr(15,M-4)): bad1+=1
    if mn < Fr(1,16)-Fr(2,M): bad2+=1
print("violations of (1/16)(1-15/(M-4)) for 9<=M<=3000:",bad1,"; of 1/16-2/M for 5<=M<=3000:",bad2)
for K in (6,7):
    for M in (50,500,3000):
        mn=min(Fr(comb(k,K)+comb(M-k,K),comb(M,K)) for k in range(M+1))
        print("k=%d M=%d min p = %.6f vs 2^(1-k)=%.6f"%(K,M,float(mn),2**(1-K)))
