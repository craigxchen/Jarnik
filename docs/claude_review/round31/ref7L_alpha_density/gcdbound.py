import sympy as sp
from sympy import I, Poly, symbols, expand, Rational as Q
from fractions import Fraction as F
from spec import tuple_at, ggcd, nrm
t=symbols('t')
def cleared_polys(n,seeds):
    A=sp.Matrix([[t**2+1,t],[t,1]])**n
    G=[]
    for s in seeds:
        a,b=s.numerator,s.denominator
        v=A*sp.Matrix([b,a]); G.append(expand(v[0]+I*v[1]))
    Gb=[expand(g.subs(I,-I)) for g in G]
    beta=t+2*I; z0=beta*sp.prod(Gb); zs=[expand(z0)]
    for (j,k) in [(0,1),(0,2),(1,2),(0,3)]:
        zs.append(expand((t-2*I)*G[j]*G[k]*sp.prod([Gb[l] for l in range(4) if l not in (j,k)])))
    return zs
for n,sd in [(1,(F(1,2),F(-2),F(7,3),F(4))),(1,(F(-1,3),F(5,2),F(6),F(-7)))]:
    zs=cleared_polys(n,sd)
    r=sp.resultant(zs[0],zs[1],t); r=sp.expand(r)
    Nr=sp.expand(r*r.subs(I,-I)); Nr=int(Nr)
    mx=0; ok=True
    for T in range(-300,301):
        z=tuple_at(n,sd,T)
        g=z[0]
        for w in z[1:]: g=ggcd(g,w)
        if g==(0,0): continue
        m=nrm(g); mx=max(mx,m); ok = ok and (Nr % m ==0)
    print(n,[str(x) for x in sd],'max N(gcd) over |T|<=300:',mx,' all divide N(Res(z0,z12)):',ok, ' N(Res) =',sp.factorint(Nr))
