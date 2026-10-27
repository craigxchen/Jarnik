# Bounds on r_0 = dim Poly_{<m'}|_A at regime sdc_q (natural logs, mpmath 50 digits).
# lower: r_0 >= #{k : |supp k| <= j} = sum_{s<=j} binom(q,s) 2^s   (functions on the Hamming ball B_j subset A are
#        polynomials of degree <= 2j < m', and restriction Poly_{<m'}|_A -> R^{B_j} is onto)
# upper: r_0 <= #{alpha in {0,1,2}^q : |alpha| < m'} = sum_{a+2b<m'} q!/(a! b! (q-a-b)!)
from mpmath import mp, mpf, log, binomial, factorial, loggamma, exp, fsum
mp.dps=50
from regime_budget import params
for q in [3251]:
    d,r,n,w,j,m=params(q)
    low=fsum(binomial(q,s)*mpf(2)**s for s in range(j+1))
    terms=[]
    for b in range(0,(m-1)//2+1):
        for a in range(0,m-2*b):
            if a+b<=q: terms.append(loggamma(q+1)-loggamma(a+1)-loggamma(b+1)-loggamma(q-a-b+1))
    M=max(terms); up=M+log(fsum(exp(t-M) for t in terms))
    U_low=q*log(3)-log(2)
    print(f"q={q} j={j} m'={m}: ln r0 in [{float(log(low)):.1f}, {float(up):.1f}] ; q ln3 = {float(q*log(3)):.1f}; "
          f"BV codim(V_H in V) < 1 needs ln H > r0 ln U/2 >= exp({float(log(low)):.1f}) ; record ln H' = m'W/40 ~ {float(m*q*296/20):.3g} (tau=296)")
