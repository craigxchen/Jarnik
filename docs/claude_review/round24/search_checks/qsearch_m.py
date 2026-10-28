# exact search: rational m = p/q (|p|,|q| <= H) where the b-quartic R(m,Z) has rational roots; check A != 0 at branch points
import sympy as sp, math, sys
m,Z=sp.symbols('m Z')
A=m**5-9*m**4+18*m**3+54*m**2+13*m-13
B=-4*m**5+36*m**4-104*m**3+72*m**2+108*m-108
D=-m**5+9*m**4-18*m**3-54*m**2+243*m-243
cub=m**3-11*m**2+39*m-13; qua=m**4+70*m**2-48*m+9
print('Res(A,cubic) =',sp.resultant(A,cub,m),' Res(A,quartic) =',sp.resultant(A,qua,m),' A(-3) =',A.subs(m,-3))
print('cubic, quartic irreducible over Q:',sp.Poly(cub,m).is_irreducible,sp.Poly(qua,m).is_irreducible)
H=int(sys.argv[1]); one=[]; full=[]
for q in range(1,H+1):
    for p in range(-H,H+1):
        if math.gcd(p,q)!=1: continue
        mv=sp.Rational(p,q)
        if mv in (1,3,-1,-3): continue
        a,b,d=[int(x*q**5) for x in (A.subs(m,mv),B.subs(m,mv),D.subs(m,mv))]
        poly=sp.Poly(a*Z**4+b*Z**2-b*Z+d,Z)
        if poly.is_zero: continue
        fl=sp.factor_list(poly)[1]
        lin=sum(e for f,e in fl if sp.Poly(f,Z).degree()==1)
        if lin>=1: one.append((mv,lin))
        if lin==4: full.append(mv)
print('H =',H,': m with >=1 rational root:',one[:30],' (count',len(one),')')
print('m with complete splitting:',full)
