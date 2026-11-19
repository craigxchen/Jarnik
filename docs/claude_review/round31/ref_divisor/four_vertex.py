from sympy import factorint
from math import isqrt, gcd
N=131490079981; ds=[6451601,28879877,92794693]
print(factorint(N)); 
for d in ds:
    s=isqrt(d); t=isqrt(d-s*s); print(d,factorint(d),s,t,t*t==d-s*s,gcd(s,t),256*t**4*N<=d*d, N%d==0)
sup=[set(factorint(d)) for d in ds]
print("sym diff", sup[2]==sup[0]^sup[1], "union=supp N", sup[0]|sup[1]==set(factorint(N)))
A=[(isqrt(d),isqrt(d-isqrt(d)**2)) for d in ds]
def mul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
c=lambda a:(a[0],-a[1])
print("Im(A1 A2 conj A3)=",mul(mul(A[0],A[1]),c(A[2])))
