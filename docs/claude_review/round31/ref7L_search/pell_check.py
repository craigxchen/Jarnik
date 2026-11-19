# independent check of Astra's polynomial Pell five-mark family (VII.5)
import sympy as sp
t=sp.symbols('t'); I=sp.I
def conjp(p): return sp.expand(sp.Poly(p,t).as_expr().subs(I,-I))
def family(n,seeds,edges):
    A=sp.Matrix([[t**2+1,t],[t,1]])
    An=A**n
    G=[]
    for s in seeds:
        v=An*sp.Matrix([1,s]); G.append(sp.expand(v[0]+I*v[1]))
    Gb=[conjp(g) for g in G]
    beta=t+2*I; betab=t-2*I
    r=len(seeds)
    z0=sp.expand(beta*sp.prod(Gb))
    Z=[z0]
    for (j,k) in edges:
        Z.append(sp.expand(betab*G[j]*G[k]*sp.prod([Gb[l] for l in range(r) if l not in (j,k)])))
    return Z
def analyse(Z):
    D=sp.degree(Z[0],t)
    N=[sp.expand(z*conjp(z)) for z in Z]
    normok=all(sp.expand(x-N[0])==0 for x in N)
    md=max(sp.degree(sp.expand(Z[i]-Z[j]),t) for i in range(len(Z)) for j in range(i))
    g=Z[0]
    for z in Z[1:]: g=sp.gcd(g,z)
    lc=[sp.Poly(z,t).LC() for z in Z]
    return dict(D=D,h=D-md,norm=normok,gcd_deg=sp.degree(g,t),distinct=len(set(map(str,Z)))==len(Z),lc=set(lc))
for n in [1,2]:
    print("5 marks n=",n,analyse(family(n,[sp.Rational(1,3),sp.Rational(2,7),sp.Rational(-5,4),sp.Rational(3,11)],[(0,1),(0,2),(1,2),(0,3)])))
    print("4 marks n=",n,analyse(family(n,[sp.Rational(1,3),sp.Rational(2,7),sp.Rational(-5,4)],[(0,1),(0,2),(1,2)])))
