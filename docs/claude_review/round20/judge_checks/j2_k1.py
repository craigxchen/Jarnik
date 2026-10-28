# Judge check J2 (exact integers): wild-killer Prop K.1. z=((b-i)^2, b^2+1, (b+i)^2), b even.
# Checks: equal norm R^2 with R=b^2+1; Gaussian gcd 1; diameter/sqrt(R) -> 4 (exact via integer chord^2);
# Pluecker coords Im(conj(z_a) z_b), their gcd = 2b; Sigma_min = log R exactly via pair-ratio heights.
from math import gcd, isqrt
from fractions import Fraction
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def nrm(a): return a[0]*a[0]+a[1]*a[1]
def gdivmod(a,b):
    n=nrm(b); num=gmul(a,(b[0],-b[1]))
    q=((2*num[0]+n)//(2*n),(2*num[1]+n)//(2*n)); r=(a[0]-gmul(q,b)[0],a[1]-gmul(q,b)[1]); return q,r
def ggcd(a,b):
    while b!=(0,0): a,b=b,gdivmod(a,b)[1]
    return a
for b in [10,100,10**4,10**8,10**20]:
    z=[gmul((b,-1),(b,-1)),(b*b+1,0),gmul((b,1),(b,1))]; R=b*b+1
    assert all(nrm(w)==R*R for w in z)
    g=ggcd(ggcd(z[0],z[1]),z[2]); assert nrm(g)==1
    d2=max(nrm((z[i][0]-z[j][0],z[i][1]-z[j][1])) for i in range(3) for j in range(i+1,3))   # = diam^2 (chord)
    ratio2=Fraction(d2,R)    # (diam/sqrt R)^2
    pl=[z[i][0]*z[j][1]-z[i][1]*z[j][0] for i in range(3) for j in range(i+1,3)]
    G=0
    for p in pl: G=gcd(G,abs(p))
    hl_num=max(abs(p) for p in pl)//G     # H(ell) = max|p|/gcd  (archimedean only after dividing by gcd)
    # pair-ratio heights h(z_x/z_a) = log(|z_x|/|gcd|) for equal-norm points: = 0.5*log(Norm(z_x)/Norm(gcd))
    def hpair(i,j): return Fraction(R*R, nrm(ggcd(z[i],z[j])))   # exp(2h)
    anchor=[ [hpair(x,a) for x in range(3) if x!=a] for a in range(3)]
    prods=[a0*a1 for a0,a1 in anchor]   # exp(2*Sigma_a)
    print(f"b={b}: (diam/sqrtR)^2={float(ratio2):.8f} (->16), Pluecker gcd={G}, 2b={2*b}, H(ell)={hl_num}, 2(b^2-1)={2*(b*b-1)}, "
          f"exp(2 Sigma_min)==R^2: {min(prods)==R*R}")
