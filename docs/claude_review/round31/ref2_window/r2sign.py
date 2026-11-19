import sympy as sp, itertools
Y=sp.symbols('Y'); ks=sp.symbols('k1:7')
F=sp.prod([Y+k for k in ks]); 
e=[1]+[sum(sp.prod(c) for c in itertools.combinations(ks,j)) for j in range(1,7)]
Q=Y**3+e[1]/2*Y**2+(e[2]/2-e[1]**2/8)*Y+e[3]/2-e[1]*e[2]/4+e[1]**3/16
R=sp.Poly(sp.expand(F-Q**2),Y)
r2=R.coeff_monomial(Y**2)
# express in centred variables: shift by mean
m=sum(ks)/6; d=sp.symbols('d1:7')
r2c=sp.factor(sp.expand(r2.subs({ks[i]:d[i] for i in range(6)}).subs(d[5],-sum(d[:5]))))
print('r2 (centred, d6=-sum):', r2c)
import random
vals=[float(r2.subs({k:random.randint(0,20) for k in ks})) for _ in range(2000)]
print('min,max r2 on random integer k:',min(vals),max(vals))
