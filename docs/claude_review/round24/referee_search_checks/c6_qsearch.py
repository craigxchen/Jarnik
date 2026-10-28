import flint, sys, math, sympy as sp
H=int(sys.argv[1])
def Gpoly(p,q):
    m=flint.fmpq(p,q)
    c=m**3-11*m**2+39*m-13
    if (m-3)==0 or (m-1)==0: return None
    u=-(m+3)*c/(16*(m-3)**2*(m-1)**2); v=-(m+1)*c/(4*(m-3)**3*(m-1))
    # (Y^2-Y+u)^2+vY = Y^4 -2Y^3 + (1+2u)Y^2 + (v-2u)Y + u^2
    return flint.fmpq_poly([u*u, v-2*u, 1+2*u, -2, 1])
hits=[]; cnt=0
for q in range(1,H+1):
    for p in range(-H,H+1):
        if math.gcd(p,q)!=1: continue
        G=Gpoly(p,q)
        if G is None: continue
        cnt+=1
        if G.degree()<4: continue
        fac=G.factor()[1]
        lin=[f for f,e in fac if f.degree()==1]
        if lin: hits.append((p,q,[str(f) for f,e in fac]))
print('H',H,'m tested',cnt,'with a rational root:',len(hits))
for h in hits[:40]: print(h)
# Galois groups at a few m
x=sp.symbols('x')
for mm in [sp.Rational(-2),sp.Rational(-3,2),sp.Rational(5),sp.Rational(7,3)]:
    c=mm**3-11*mm**2+39*mm-13; u=-(mm+3)*c/(16*(mm-3)**2*(mm-1)**2); v=-(mm+1)*c/(4*(mm-3)**3*(mm-1))
    G=sp.Poly(sp.expand((x**2-x+u)**2+v*x),x)
    print('m=',mm,' Galois group over Q:',sp.polys.numberfields.galoisgroups.galois_group(G,by_name=True)[0])
