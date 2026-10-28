# Judge check J2: for nested blocks, gcd(A(x),C(x)) can exceed prod d_T^{mu_T} by a full block.
# Unit pencil A=(x1-x2)(x4-x5), C=(x1-x3)(x2-x3) (anchored chart, x0=infinity). Chain at p: {1,2} subset {1,2,3}.
from math import gcd
from itertools import combinations
import sympy
p=10007
x=[None,0,p*p*3,p*5,11,13*17]   # x1=0, x2 = 3p^2 (so p^2 | x1-x2), x3 = 5p (p || x1-x3, x2-x3)
A=(x[1]-x[2])*(x[4]-x[5]); C=(x[1]-x[3])*(x[2]-x[3])
def vp(n,p):
    n=abs(n);v=0
    while n%p==0: n//=p; v+=1
    return v
# contacts at p: d_{123} gets p^1 (all of 1,2,3 agree mod p), d_{12} gets p^1 (1,2 agree mod p^2 beyond)
nu={'12':{'A':1,'C':0},'123':{'A':1,'C':2}}
mu={T:min(v['A'],v['C']) for T,v in nu.items()}
print("v_p(A)=",vp(A,p)," v_p(C)=",vp(C,p)," v_p(gcd)=",vp(gcd(A,C),p))
print("sum_T mu_T v_p(d_T) =", mu['12']*1+mu['123']*1, " -> excess", vp(gcd(A,C),p)-(mu['12']+mu['123']), "copies of p = one full block d_12")
