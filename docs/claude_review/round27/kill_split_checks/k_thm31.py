# Independent check of split.md Theorem 3.1 (mnn31 curve, quadric cone). Own code.
import sympy as sp, random, itertools, math
from fractions import Fraction as Fr
a,b,c,d,t,u,w=sp.symbols('a b c d t u w')
k=a*b/(c*d)
U=(k*t*(c+d)-(a+b))/(1-k*t**2); Wt=t*U
conic=sp.simplify(U*(U+a+b)-k*Wt*(Wt+c+d))
print("A conic residual:",conic)
D=[(u+a)*(u+b),(u+a)*(u+a+b),(u+b)*(u+a+b),k*(w+c)*(w+c+d),k*(w+d)*(w+c+d)]
Dt=[sp.together(Di.subs({u:U,w:Wt})) for Di in D]
print("A (u+a)(u+b)-k(w+c)(w+d) on conic:",sp.simplify(Dt[0]-(k*(w+c)*(w+d)).subs({u:U,w:Wt})))
Y=[sp.together(Dt[j]/Dt[0]-1) for j in range(5)]
print("A cone identity Y2Y3-Y4Y5:",sp.simplify(Y[1]*Y[2]-Y[3]*Y[4]))
# B: valuation vectors at all places, for several parameter values
def places(vals):
    sub=dict(zip((a,b,c,d),vals))
    Ds=[sp.factor(sp.together(Di.subs(sub))) for Di in Dt]
    facs=set()
    for Di in Ds:
        n_,d_=sp.fraction(sp.together(Di))
        for f,_ in sp.factor_list(n_,t)[1]+sp.factor_list(d_,t)[1]:
            if f.has(t): facs.add(sp.expand(f))
    res=[]
    for f in facs:
        v=[]
        for Di in Ds:
            n_,d_=sp.fraction(sp.together(Di))
            on=0
            while sp.rem(n_,f,t)==0: n_=sp.quo(n_,f,t); on+=1
            od=0
            while sp.rem(d_,f,t)==0: d_=sp.quo(d_,f,t); od+=1
            v.append(on-od)
        m=min(v); v=[x-m for x in v]
        res.append((sp.degree(f,t),str(f),tuple(v)))
    v=[]
    for Di in Ds:
        n_,d_=sp.fraction(sp.together(Di)); v.append(sp.degree(d_,t)-sp.degree(n_,t))
    m=min(v); res.append((1,'inf',tuple(x-m for x in v)))
    return res
for vals in [(2,3,1,2),(sp.Rational(3,7),sp.Rational(5,2),sp.Rational(2,3),sp.Rational(11,5))]:
    print("B params",vals,"k=",sp.Rational(vals[0]*vals[1])/(vals[2]*vals[3]))
    cuts=[]; BC=0
    for dg,f,v in sorted(places(vals)):
        if any(v):
            B=frozenset(i+1 for i in range(5) if v[i]>0)
            assert set(v)<={0,1}, v
            cuts.append(B); BC+=dg
        print("   place deg",dg,f,v)
    canon=[min(B,frozenset(range(1,6))-B,key=lambda s:(len(s),sorted(s))) for B in cuts]
    print("   B.C =",BC," distinct cuts:",len(set(canon)),"of",len(canon))
# C: Jacobian rank of (a,b,c,d,t)->(Y2..Y5) exact
F=sp.Matrix([Y[j] for j in range(1,5)]); J=F.jacobian([a,b,c,d,t])
random.seed(1)
for trial in range(4):
    pt={s:sp.Rational(random.randint(1,40),random.randint(1,13)) for s in (a,b,c,d,t)}
    print("C Jacobian rank at",list(pt.values()),":",J.subs(pt).rank())
