from math import factorial, comb
import mpmath as mp
mp.mp.dps=50
# Lemma 5.6' at j=91: (9j+72)^63 * 2^129 < (2j-30)^63 * 63!
for j in (90,91,92,200):
    lhs=(9*j+72)**63*2**129; rhs=(2*j-30)**63*factorial(63)
    print(j, lhs<rhs)
# sharper: binom(e+63,63) vs (2j-30)^63/2^129 at e=9j+9
for j in (60,70,80,91):
    print("binom check j",j, comb(9*j+9+63,63) < (2*j-30)**63//2**129)
# I_t and asymptotic rho/j via Irwin-Hall: need P(sum_{b} U[0,1] <= x/2) >= 1/I_t, x = rho/j
def IH_cdf(b,s):
    s=mp.mpf(s)
    return sum((-1)**k*mp.binomial(b,k)*(s-k)**b for k in range(0,int(mp.floor(s))+1))/mp.factorial(b)
for t in range(2,8):
    b=2**t-1; It=mp.mpf(2)**(t-1+(t-2)*(2**(t-1)-1))
    lo,hi=mp.mpf(0),mp.mpf(b)
    for _ in range(200):
        mid=(lo+hi)/2
        if IH_cdf(b,mid)>=1/It: hi=mid
        else: lo=mid
    print("t",t,"b",b,"log2 I_t",t-1+(t-2)*(2**(t-1)-1),"rho/j >= ",mp.nstr(2*hi,6))
