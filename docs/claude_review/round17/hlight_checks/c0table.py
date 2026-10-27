# Rigorous lower bound LB(q) for sigma over all non-pair integral lambda (hlight.md Thm 2.1), exact rationals;
# sqrt(q) bounded above by an integer-rounded rational (isqrt+1). Prints c0(q)=LB/(2q).
from fractions import Fraction as Fr
from math import isqrt
import sympy, sys
def g(n,q): return Fr(q*(n-2)+4*n-n*n,n+2)
def LB(q):
    sq=Fr(isqrt(q)+1)                      # >= sqrt(q)
    b4=max(Fr(q,3), Fr(q,2)-11-Fr(3,2)*sq) # F=0, n=4
    bF1=Fr(q-1,2)                          # |F|=1 (n=4 sharp; n>=6 gives (q+1)/2)
    bF0=min([max(Fr(n-1),g(n,q)) for n in range(6,(q+1)//2+1,2)]+[Fr(10**9)])  # F=0, 6<=n<=(q+1)/2
    bF0big=Fr((q+1)//2+1-1) if True else None  # F=0, n>(q+1)/2, n even => n>=(q+3)/2, sigma>=n-1>=(q+1)/2
    bnt=Fr(q-1)                            # non-ternary or |F|>=2
    parts={'F0n4':b4,'F1':bF1,'F0n6+':bF0,'F0big':Fr(q+1,2),'nonternary':bnt}
    k=min(parts,key=parts.get)
    return parts[k],k,parts
for q in [int(a) for a in sys.argv[1:]] or [7,11,19,23,31,43,83,199,499,1019,2003,3251,5479,10007,100003]:
    if not sympy.isprime(q) or q%4!=3: continue
    lb,k,p=LB(q)
    import math
    print(f"q={q} LB={float(lb):.2f} (ceil {math.ceil(lb)}) c0=LB/r={float(lb/(2*q)):.4f} binding={k}  [F0n4={float(p['F0n4']):.1f}, F1={float(p['F1']):.1f}, F0n6+={float(p['F0n6+']):.1f}]",flush=True)
