# [J] judge checks for Part IV: ledger algebra, thresholds, trace seed (GF(2^d))
import sympy as sp
from fractions import Fraction as Fr
m,s,L,d,E,S=sp.symbols('m s L delta E S')
# per-layer: S coefficient 8/9+(4/9)(m-1)/2, S <= sum s_h - mL/4 + m delta
coef=sp.Rational(8,9)+sp.Rational(4,9)*(m-1)/2
print('S coeff', sp.factor(coef))
perlayer=coef*(s-m/4)-sp.Rational(4,9)*s*(s-1)/2+4-sp.Rational(8,9)*s
print('per-layer identity ok:', sp.simplify(perlayer-(4-m/6-sp.Rational(2,9)*(s-m/2)**2))==0)
tot=sp.Rational(9,2)*(L*(4-m/6))
print('9/2*(4-m/6) =', sp.factor(sp.Rational(9,2)*(4-m/6)), ' ; target -(3/4)(m-24):', sp.expand(-sp.Rational(3,4)*(m-24)))
print('delta coeff*9/2:', sp.expand(sp.Rational(9,2)*coef*m), ' E coeff:', sp.Rational(9,2)*sp.Rational(8,9))
# phi superadditivity
phi=lambda k:max(0,k*(k-1)//2+2*k-9)
print('phi defect ok:', all(phi(u+v)-phi(u)-phi(v)<=2*u*v for u in range(1,80) for v in range(1,80)))
# sample-size constants
print('48 marks:', Fr(3*24,4*48*51), ' 49 marks:', Fr(3*25,4*49*52))
best=max(((Fr(3*(k-24),4*k*(k+3)),k) for k in range(25,200)))
print('best sampled excess', best, float(best[0]))
# FF thresholds: L=2D, delta<=c ; m marks impossible iff 3(m-24)D > 2m(m+3)c
for mm in (25,48,49):
    print('m=%d excluded for D > %s c (=%.3f c)'%(mm,Fr(2*mm*(mm+3),3*(mm-24)),2*mm*(mm+3)/(3*(mm-24))))
# odd D, c=0: delta=-1/2: m=24 gives sum sq <= -m(m+3)/2 <0
print('odd D, m=24 RHS', -24*27/2)
# trace seed: LHS of (IV.2) with delta=E=0, m=L=q-1, s=q/2-1
for q in (8,16,32,64,128):
    mq=q-1; lhs=Fr(3,4)*(mq-24)*mq+mq*(Fr(q,2)-1-Fr(mq,2))**2
    print('q=%d LHS=%s'%(q,lhs))
# independent GF(2^d) check of the trace seed degrees and divisibility
def gfmul(a,b,d,poly):
    r=0
    while b:
        if b&1:r^=a
        b>>=1;a<<=1
        if a>>d:a^=poly
    return r
def gfpow(a,e,d,poly):
    r=1
    while e:
        if e&1:r=gfmul(r,a,d,poly)
        a=gfmul(a,a,d,poly);e>>=1
    return r
def tr(x,d,poly):
    t=0;y=x
    for _ in range(d):t^=y;y=gfmul(y,y,d,poly)
    return t
for d_,poly in ((3,0b1011),(4,0b10011),(5,0b100101),(6,0b1000011)):
    q=1<<d_; inv=lambda a:gfpow(a,q-2,d_,poly)
    def fa(a): # coefficients of f_a: t^(2^j-1) coeff a^(2^j - q/2), j<=d-2
        c=[0]*(q//2)
        for j in range(d_-1):
            c[(1<<j)-1]=gfpow(a,((1<<j)-q//2)%(q-1),d_,poly)
        return c
    def ev(c,x):
        r=0;p=1
        for co in c:r^=gfmul(co,p,d_,poly);p=gfmul(p,x,d_,poly)
        return r
    ok=True
    for a in range(1,q):
        Ma=[x for x in range(1,q) if tr(gfmul(a,x,d_,poly),d_,poly)==0]
        assert all(ev(fa(a),x)==gfpow(x,q//2-1,d_,poly) for x in Ma)
        for b in range(a+1,q):
            diff=[u^v for u,v in zip(fa(a),fa(b))]
            deg=max(i for i,v in enumerate(diff) if v)
            common=[x for x in Ma if tr(gfmul(b,x,d_,poly),d_,poly)==0]
            ok&= deg==q//4-1==len(common) and all(ev(diff,x)==0 for x in common)
    print('q=%d trace seed: deg(f_a-f_b)=q/4-1=#common roots, all vanish:'%q, ok)
