import sympy as sp
from fractions import Fraction as Fr
import sys; sys.path.insert(0,'.')
from ref_beta7 import setup
lam=sp.symbols('lam')
T9=[(2,3,-5,1),(Fr(1,3),3,2,Fr(5,3)),(-2,Fr(7,5),-3,1),(2,Fr(9,5),Fr(6,5),1),(Fr(5,2),Fr(-9,5),-7,1),
       (Fr(-3,4),Fr(-1,3),Fr(1,2),1),(-3,Fr(1,2),4,1),(Fr(3,5),2,Fr(-3,2),1),(Fr(7,5),Fr(3,2),Fr(7,3),1)]
for t in T9:
    S=setup(t); e1,e2,e3,e4=[sp.Rational(str(v)) for v in S['e']]
    s=sp.Rational(str(S['s'])); tp=sp.Rational(str(S['tp']))
    def Dl(TT):
        a2=1+lam*e2+lam*TT; a1=-lam*(e3+e1*TT); a0=lam*e4-TT
        return sp.expand(a1**2-4*a2*a0)
    ds=sp.Poly(Dl(s),lam); dt=sp.Poly(Dl(tp),lam)
    r=sp.resultant(ds,dt)
    lb=sp.Rational(str(S['lb'])); lc=sp.Rational(str(S['lc']))
    rs=sp.roots(ds); rt=sp.roots(dt)
    # degenerate members: common root of P_lam, Q_lam in z
    z=sp.symbols('z')
    P=(1+lam*e2)*z**2-lam*e3*z+lam*e4; Q=1+lam*e1*z-lam*z**2
    rdeg=sp.factor(sp.resultant(P,Q,z))
    print([str(v) for v in t],'Res(Delta(s),Delta(t\'))!=0:',r!=0,'roots Delta(s):',list(rs),'lam_b',lb,'roots Delta(t\'):',list(rt),'lam_c',lc,'| deg-res:',rdeg)
