import flint, sympy, pickle, itertools
from param import load
from ms import system
from cuttypes import CUTS
from gb import cuts_to_columns
par=load('p.out'); p=par['p']; w=par['w']
T=pickle.load(open('types_8_4.pkl','rb'))
v,e=system(cuts_to_columns(T[0],CUTS),4,2,1)
syms=sympy.symbols(v)
polys=[sympy.Poly(sympy.sympify(eq.replace('^','**')),*syms) for eq in e]
def evalp(P,X):
    val=flint.nmod_poly([0],p)
    for mon,c in P.terms():
        t=flint.nmod_poly([int(c)%p],p)
        for nm,ex in zip(v,mon):
            if ex: t=(t*pow(X[nm],ex))%w
        val=(val+t)%w
    return val.degree()<0
for order in ['fwd','rev']:
  nums=par['nums'] if order=='fwd' else par['nums'][::-1]
  for sgn in [1,-1]:
    for deriv in [True,False]:
      dw=w.derivative() if deriv else flint.nmod_poly([1],p)
      g,a,b=w.xgcd(dw); inv=(b*pow(int(g.coeffs()[0]),p-2,p))%w
      X={nm:(sgn*num*inv)%w for nm,num in zip(v[:-1],nums)}
      X['u']=flint.nmod_poly([0,1],p)
      ok=[evalp(P,X) for P in polys]
      print(order,sgn,deriv,ok)
