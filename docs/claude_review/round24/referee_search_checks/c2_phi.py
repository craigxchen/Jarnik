import sympy as sp, pickle, flint
b0,b1,b2,e2,e3,e4,u,v,m,k,W = sp.symbols('b0 b1 b2 e2 e3 e4 u v m k W')
facs=pickle.load(open('c1b_factors.pkl','rb'))
F24=[f for f,e in facs if e==1 and len(f)>500][0]
F24=sp.sympify(F24.replace('^','**'))
b3=-(b0+b1+b2); bs=[b0,b1,b2,b3]
E2=sp.expand(sum(bs[i]*bs[j] for i in range(4) for j in range(i+1,4)))
E3=sp.expand(sum(bs[i]*bs[j]*bs[k_] for i in range(4) for j in range(i+1,4) for k_ in range(j+1,4)))
E4=sp.expand(b0*b1*b2*b3)
Fclaim=32*e2**6*e4**3-96*e2**5*e3**2*e4**2+37*e2**4*e3**4*e4-256*e2**4*e4**4-4*e2**3*e3**6+1024*e2**3*e3**2*e4**3-992*e2**2*e3**4*e4**2+288*e2*e3**6*e4-27*e3**8-256*e3**4*e4**3
Fb=sp.expand(Fclaim.subs({e2:E2,e3:E3,e4:E4}))
q=sp.cancel(F24/Fb)
print('F24 / F(e2,e3,e4) =',q)
# Phi(u,v): F/e2^12 with e4=u e2^2, e3^2 = v e2^3
Phi=sp.expand(sp.expand(Fclaim.subs(e4,u*e2**2)).subs(e3,sp.sqrt(v*e2**3))/e2**12)
Phi=sp.factor(Phi); print('Phi =',Phi)
c=m**3-11*m**2+39*m-13
um=-(m+3)*c/(16*(m-3)**2*(m-1)**2); vm=-(m+1)*c/(4*(m-3)**3*(m-1))
print('Phi(u(m),v(m)) =',sp.simplify(Phi.subs({u:um,v:vm})))
print('Phi irreducible over Q:',sp.factor_list(sp.Poly(Phi,u,v))[1].__len__(), 'factors')
# triple point & pencil
Pk=sp.factor(sp.expand(Phi.subs(v,k*u)))
print('Phi(u,ku) =',Pk)
quad=sp.cancel(Pk/u**3); qp=sp.Poly(quad,u)
print('deg in u',qp.degree(),' disc =',sp.factor(qp.discriminant()))
# birationality of m -> (u,v): degree of the map = number of m' with (u(m'),v(m')) = (u(m),v(m)) generically
m2=sp.symbols('m2')
g1=sp.numer(sp.together(um-um.subs(m,m2))); g2=sp.numer(sp.together(vm-vm.subs(m,m2)))
G=sp.gcd(sp.Poly(g1,m2),sp.Poly(g2,m2))
print('gcd in m2:',sp.factor(G.as_expr()))
# irreducible over C? check absolute irreducibility via factoring over Q(sqrt?) -- genus 0 rational param implies irreducible image anyway
pickle.dump((str(Phi),),open('c2_phi.pkl','wb'))
