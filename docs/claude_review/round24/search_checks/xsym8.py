import sympy as sp, pickle, numpy as np, itertools
u,v,k,m,W=sp.symbols('u v k m W')
Phi=-256*u**4 - 256*u**3*v**2 + 1024*u**3*v + 32*u**3 - 992*u**2*v**2 - 96*u**2*v + 288*u*v**3 + 37*u*v**2 - 27*v**4 - 4*v**3
G=sp.Poly(sp.expand(Phi.subs(v,k*u)/u**3),u)
c2,c1,c0=[G.coeff_monomial(u**d) for d in (2,1,0)]
kk=4*(m**2-1)/(m**2-9); sq=sp.factor((kk-4)**2*(9*kk-4)*32*m/(m**2-9))   # sqrt(Delta) = (k-4)^2 (9k-4) sqrt((k-4)(9k-4))
for sgn in (1,-1):
    U=sp.factor(sp.simplify(((-c1+sgn*sq)/(2*c2)).subs(k,kk)))
    V=sp.factor(sp.simplify(kk*U))
    print('branch',sgn,': u(m) =',U,'\n            v(m) =',V)
    print('   check Phi(u(m),v(m)) =',sp.simplify(Phi.subs({u:U,v:V})))
    pickle.dump((U,V),open(f'uv_branch{sgn}.pkl','wb'))
# numeric family points -> (u,v)
import xstruct as X
P=np.load('xfam_points.npy')
for p in P[[0,150,300]]:
    rho,r4,sig,beta,q=X.unpack(p)
    bs=[1.0,rho[0].imag,rho[1].imag,rho[2].imag]
    e2=sum(bs[i]*bs[j] for i,j in itertools.combinations(range(4),2)); e3=sum(bs[i]*bs[j]*bs[l] for i,j,l in itertools.combinations(range(4),3)); e4=np.prod(bs)
    uu=e4/e2**2; vv=e3**2/e2**3; kv=vv/uu
    # m from k: k(m^2-9) = 4(m^2-1) -> m^2 (k-4) = 9k-4
    m2=(9*kv-4)/(kv-4)
    print('family point: sum b =',round(sum(bs),12),' u =',uu,' v =',vv,' Phi =',float(Phi.subs({u:uu,v:vv})),' m^2 =',m2)
