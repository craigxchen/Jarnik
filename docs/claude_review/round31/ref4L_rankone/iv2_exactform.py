# Re-read Part IV.2 models as EXACT clusters for their interpolating rational form Q_*.
from fractions import Fraction as Fr
from sympy import factorint, legendre_symbol, isprime
import itertools, math
x=[-27,-270,135,120,40]; s1,s2,s3,s5=-2,-53925,-5042250,4723920000; l=97200
b0=3*s5//2
def model(P,q):
    L=l*q; B=b0*q*q
    h=[xi**3-s1*xi**2+s2*xi for xi in x]
    U=[P+q*hi for hi in h]
    V=[(P*U[i]+B*x[i], L*U[i]) for i in range(5)]
    J={}
    for i,j in itertools.combinations(range(5),2):
        J[(i,j)]=P+q*x[i]*x[j]*(s1-x[i]-x[j])
    n=[math.prod(J[tuple(sorted((i,j)))] for j in range(5) if j!=i) for i in range(5)]
    return V,J,n
def solve3(V,n,idx):
    # Q=aX^2+bXY+cY^2 through three points
    import sympy as sp
    M=sp.Matrix([[V[i][0]**2,V[i][0]*V[i][1],V[i][1]**2] for i in idx]); r=sp.Matrix([n[i] for i in idx])
    return list(M.LUsolve(r))
import sympy as sp
for (P,q) in [(10**15+37,30),(10**25+7,30),(10**40+9,210)]:
    V,J,n=model(P,q)
    a,b,c=solve3(V,n,[0,1,2])
    ok=all(a*X*X+b*X*Y+c*Y*Y==n[i] for i,(X,Y) in enumerate(V))
    detQ=a*c-b*b/4; disc=b*b-4*a*c
    den=sp.ilcm(*[sp.fraction(t)[1] for t in (a,b,c)])
    D=sp.Rational(disc)*den**2
    Ent=[1-a,-b,1-c]
    detE=(1-a)*(1-c)-b*b/4
    print('P~10^%d q=%d exact on all five:'%(len(str(P))-1,q),ok,' den(Q*)=',den,
          ' disc(Q*) sign',sp.sign(disc),' log10|disc|=%.2f'%math.log10(abs(float(disc))),' log10 P=%.2f'%math.log10(P),
          ' det E=',detE, ' disc*den^2 is square:',sp.sqrt(D).is_rational)
    if P<10**16:
        Dint=int(D)
        bad=0;tot=0
        for e,Je in J.items():
            for p in factorint(Je):
                if p==2: continue
                tot+=1
                if Dint%p==0: print('  ramified block prime',p)
                elif legendre_symbol(Dint%p,p)!=1: bad+=1
        print('  block prime factors checked',tot,' inert:',bad)
