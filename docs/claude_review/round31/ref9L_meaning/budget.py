# Entropy budgets for the random-parent theorem (IX.0) and for fixed parents. Own code.
from math import log, lgamma, exp, ceil
def lbinom(n,k):
    if n>1e12*k: return k*log(n/k)+k-0.5*log(2*3.141592653589793*k)
    return lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1)
print("== IX.0 numbers: pool of K primes, s drawn, J=log(K/s)")
for s in [10**3,10**4,10**6,10**9]:
    for J in [40,60]:
        K=s*exp(J); P=K*log(K)  # dyadic pool [P,2P] has ~P/log P primes (all primes); split: halve
        Lpar=lbinom(K,s)            # log #parents
        Lexc=Lpar - s*J/20          # log #exceptional parents allowed
        b=1+ceil(100*log(2*P)/J)
        W=s*log(P); logR=W/2
        print(f"s={s:.0e} J={J}: P~{P:.2e} (P/(s log s)={P/(s*log(s)):.2e}); log#parents={Lpar:.3e}, "
              f"log#exc<={Lexc:.3e} ({Lexc/Lpar:.3f} of it); b-1={b-1}, log s={log(s):.1f}; "
              f"exc fraction = R^-{(s*J/20)/logR:.3f}; logR/loglogR={logR/log(logR):.3e}")
print("== union-bound exponent per slot (twin, 4 labels; circle 7 labels), with b-1=100 log2P/J")
for J in [18,18.5,20,40]:
    for lab in [4,7]:
        print(J,lab, "exponent/s =", round(log(lab)+9*J/200 - J/5,3), " need <= -J/20 =", -J/20)
print("== fixed parent, many marks: pattern entropy vs total parent entropy")
for s in [10**3,10**6]:
    for J in [0,1,40]:
        Lpar = s*(1+J) if J>0 else 0.0
        for m in [3,10,int(log(s)),int(s**0.5)]:
            patt = s*m*log(2)        # sign patterns of m marks on s primes
            colms = lbinom(s+2**(m-1),s) if m<60 else float('inf')  # column multisets
            print(f"s={s:.0e} J={J} m={m}: log#parents={Lpar:.2e}, log#patterns={patt:.2e}, log#column-multisets={colms:.2e}")
