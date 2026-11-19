# For the Part V polynomial families: does an interpolating quadratic Q_* with Q_*(V_i)=n_i exist (IV.5)?
# Compute E = I - Q_* through marks 1,2,3 and test at 4,5.  Then the integer 191-digit instance.
import sympy as sp
from itertools import combinations
from fractions import Fraction
exec(open('defects.py').read().split('# ---- four-star')[0])   # helpers det,nrm
T=sp.symbols('T',real=True)
def fam_star():
    v,t,k=sp.Rational(3,2),2,2
    r=t**2; f=(1-r)*(r**2-3*r+1); x=(k+sp.Rational(f)/k)/2; y=(sp.Rational(f)/k-k)/(2*t)
    w=[v*(1+x)/2, v*(1-x)/2, v*(-r+y)/2, v*(-r-y)/2]
    s1=sum(w); G=sp.expand(sp.prod([T-sp.I*a for a in w])); A=-T/s1+sp.I
    u=[1]+[a**2 for a in w]
    F=[sp.expand(G+sp.Rational(u[i]-1)/u[i]*(u[i]+T**2)*A) for i in range(5)]
    J={(i,j):T**2+u[i]*u[j] for i,j in combinations(range(5),2)}
    return F,J
def test(F,J,name):
    V=[(sp.re(f).expand(),sp.im(f).expand()) for f in F]
    n=[sp.expand(sp.prod([J[e] for e in J if i in e])) for i in range(5)]
    a,b,c=sp.symbols('a b c')
    Q=lambda v: a*v[0]**2+b*v[0]*v[1]+c*v[1]**2
    sol=sp.solve([Q(V[i])-n[i] for i in range(3)],[a,b,c],dict=True)[0]
    res=[sp.factor(sp.together(Q(V[i]).subs(sol)-n[i])) for i in (3,4)]
    print(name,': Q_* through marks 1-3; residual at marks 4,5 =',[sp.Poly(sp.numer(rr),T).degree() if rr!=0 else 0 for rr in res],'(0 means interpolation identity holds)')
    print('  coefficient degrees of Q_* (num/den):',[(sp.Poly(sp.numer(sp.together(sol[s])),T).degree(),sp.Poly(sp.denom(sp.together(sol[s])),T).degree()) for s in (a,b,c)])
test(*fam_star(),'four-star')
