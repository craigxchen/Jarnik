# integer polynomial template from the exact cube: roots a_k = m*r_k in Z[i]; z_j(u) in Z[i][u]; arc constant C at integer u.
import sympy as sp, itertools
from fractions import Fraction as Fr
from math import gcd
y=[Fr(-9,4),Fr(-2,3),Fr(1,4),Fr(8,3)]; eps=[1,-1,-1,1]
r=[(Fr(e)/v, v) for e,v in zip(eps,y)]
m=1
for (x,yy) in r:
    for q in (x,yy): m=m*q.denominator//gcd(m,q.denominator)
roots=[]
for (x,yy) in r:
    a=complex(int(x*m),int(yy*m)); roots+= [a,-a]
print("scale m=",m,"Gaussian integer roots:",roots)
u=sp.symbols('u')
G=[sp.Integer(int(z.real))+sp.I*int(z.imag) for z in roots]
Q=sp.expand(sp.prod([u-g for g in G]))
cl=[]
for mask in range(256):
    z=sp.expand(sp.prod([u-sp.conjugate(g) if (mask>>k)&1 else u-g for k,g in enumerate(G)]))
    if sp.Poly(sp.expand(z-Q),u).degree()<=4: cl.append(z)
print("points:",len(cl))
import mpmath as mp
for U in (10**6,10**9,10**12):
    vals=[complex(sp.N(z.subs(u,U),40)) for z in cl]
    Zs=[sp.expand(z.subs(u,U)) for z in cl]
    # gaussian gcd of values
    g=Zs[0]
    def ggcd(a,b):
        a=complex(a);b=complex(b)
        return None
    R=abs(complex(sp.N(Zs[0],50)))
    chords=[abs(complex(sp.N(a-b,50))) for a,b in itertools.combinations(Zs,2)]
    print(U," max chord / sqrt(R) =", max(chords)/mp.sqrt(R), " (R=|z(u)|, no gcd reduction)")
