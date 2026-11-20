# Exact-rational check of the IX.4 parameter system (base-q exponents).
from fractions import Fraction as Fr
import sympy as sp

xi,ka,ta,a,ep,ga,nu,lam=sp.symbols('xi kappa tau a epsilon gamma nu lambda',real=True)
d=2-xi; b=sp.Rational(3,2)+a+ka+2*ep
# per-token exponent of primitive balanced words at (h,k): count (1-ga)(1-nu), digit density b*nu,
# atoms: alpha per erased digit, prob P^{-h/2} -> -d/2
for name,alpha in [('strong',a+sp.Rational(1,2)+ga/2),('weak',a+sp.Rational(1,2))]:
    F=(1-ga)*(1-nu)+(b-alpha)*nu-d/2
    print(name,'coef nu:',sp.simplify(sp.diff(F,nu)),' F(nu=1/2)=',sp.simplify(F.subs(nu,sp.Rational(1,2))))
# nonprimitive: count (1-ka)/2 (at h>H), digits b/4, prob -d/2
print('nonprim (at h=H):',sp.simplify((1-ka)/2+b/4-d/2))
print('unbalanced h>H_u:',sp.simplify(1-lam-d/2))
print('zeta tail per token at h>=q/2 (J<=q/4,h>=q/2):',sp.simplify(b/4-d/4 - 0) ,' (times q/h<=2 -> (b-d)/2*... must be <0)')

def feasible(X,amin=Fr(1,10**6)):
    # strict system: ka>2xi+4ep, ta>ka+xi+2ep, ta<1/2, b<d, xi/2<lam<ka/2, a>0, ep>0
    ep=Fr(1,10**7); a=amin
    ka=2*X+5*ep; ta=ka+X+3*ep
    lam=(X/2+ka/2)/2
    b=Fr(3,2)+a+ka+2*ep; d=2-X
    conds=dict(ka=ka>2*X+4*ep, ta=ta>ka+X+2*ep, tah=ta<Fr(1,2), bd=b<d, lam=(X/2<lam<ka/2), tau_vs_atom=(2*ta-2*a<1))
    Fs=X/2+0 # placeholder
    Fstrong=(X+ka+2*ep)/2-Fr(3,4)*ka   # worst gamma=kappa
    Fweak=(X+ka+2*ep-ta)/2             # worst gamma=tau
    nonp=(1-ka)/2+b/4-d/2
    return conds,dict(Fstrong=float(Fstrong),Fweak=float(Fweak),nonprim=float(nonp),d_minus_b=float(d-b),half_minus_tau=float(Fr(1,2)-ta))
for X in [Fr(1,12),Fr(1,6)-Fr(1,1000),Fr(1,6)-Fr(1,10**5),Fr(1,6)]:
    c,v=feasible(X); print(X,all(c.values()),c,v)
# LP: is the closed system feasible at xi=1/6 in the limit? sup of xi subject to ka>=2xi, ta>=ka+xi, ta<=1/2, 3/2+a+ka<=2-xi, a>=0
from scipy.optimize import linprog
# vars xi,ka,ta,a ; maximize xi
res=linprog(c=[-1,0,0,0],A_ub=[[2,-1,0,0],[1,1,-1,0],[0,0,1,0],[1,1,0,1]],b_ub=[0,0,0.5,0.5],bounds=[(0,None)]*4)
print('LP sup xi (closed system, a>=0):',-res.fun,res.x)
res=linprog(c=[-1,0,0,0],A_ub=[[2,-1,0,0],[1,1,-1,0],[0,0,1,0],[1,1,0,1],[1,1,-1+0,-1]],b_ub=[0,0,0.5,0.5,0],bounds=[(0,None)]*4)
print('LP with old a>=tau:',-res.fun,res.x)
# residue slack variant (work notes): b0<=1.01/d at d=11/6
D=Fr(11,6); print('1.01/d at 11/6 =',float(Fr(101,100)/D),' 4/3-.56/2=',float(Fr(4,3)-Fr(28,100)),' (16/9-.56)/2+4/9=',float((Fr(16,9)-Fr(56,100))/2+Fr(4,9)), ' 1-.56/2=',float(1-Fr(28,100)))
print('limits: d->',Fr(11,6),' 4d ->',4*Fr(11,6))
