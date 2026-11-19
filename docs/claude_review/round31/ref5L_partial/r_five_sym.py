# Independent check of the five-block curve (Part V.3), with own conventions [X,Y]=ReX ImY - ImX ReY.
import sympy as sp
from itertools import combinations
x=sp.symbols('x',real=True); T=sp.symbols('T',real=True)
I=sp.I
z=[sp.Rational(3,2),sp.Integer(3),sp.Rational(1,3),sp.Rational(-32,15),sp.Rational(-5,54)]
u=[a*a for a in z]; E=list(combinations(range(5),2)); h={e:z[e[0]]*z[e[1]] for e in E}
print('h values',[h[e] for e in E]); assert len({abs(h[e]) for e in E})==10
def RI(F):
    F=sp.expand(F); return sp.expand((F+sp.conjugate(F))/2), sp.expand((F-sp.conjugate(F))/(2*I))
def br(F,H):
    a,b=RI(F); c,d=RI(H); return sp.expand(a*d-b*c)
def nrm(F):
    a,b=RI(F); return sp.expand(a*a+b*b)
G=sp.expand(sp.prod([x+I*h[0,j] for j in range(1,5)]))
B=(4500480*x**2+1413120)/sp.Integer(1226641)+I*(x**2+sp.Rational(26758,283071)*x+1)
print('[G,B]-(x-1)^6 =',sp.factor(br(G,B)-(x-1)**6))
# derive d independently
c=lambda j,X: (u[j]-u[0])*(u[0]+X**2/u[j])
dd=sp.nsimplify(sp.simplify(-G.subs(x,I)/(B.subs(x,I)*c(1,I))))
print('d computed =',dd,' claimed',-sp.Rational(418284581,1009043070), ' c2(i)==c3(i):',sp.simplify(c(1,I)-c(2,I))==0)
d=-sp.Rational(418284581,1009043070); assert dd==d
Lc=[sp.expand(G+d*B*c(j,x)) for j in range(5)]
iso=set()
for (j,k) in E:
    blk=x**2+h[j,k]**2
    q=sp.cancel(br(Lc[j],Lc[k])/((x-1)**6*blk)); assert q.is_Rational and q!=0
    assert q==d*u[0]*(1/u[j]-1/u[k])
    for m in (j,k):
        if sp.rem(nrm(Lc[m]),blk,x)==0: iso.add(((j,k),m))
print('x-chart isotropic incidences:',sorted(iso))
v1=[sp.expand(Lc[j].subs(x,1)) for j in range(5)]
rat=[sp.simplify(v1[j]/v1[0]) for j in range(5)]
print('L_j(1)/L_1(1) =',rat)
F=[sp.expand(sp.cancel(sp.expand(T**4*Lc[j].subs(x,1+1/T))/v1[j])) for j in range(5)]
for f in F:
    p=sp.Poly(f,T); assert p.degree()==4 and p.LC()==1
    a,b=RI(f); assert sp.gcd(a,b)==1 or sp.degree(sp.gcd(a,b),T)==0
iso2=set(); cof={}
for (j,k) in E:
    blk=sp.expand(((T+1)**2+h[j,k]**2*T**2)/(1+h[j,k]**2))
    a_=-1/(1+h[j,k]**2); b_=h[j,k]/(1+h[j,k]**2)
    assert sp.expand((T-a_)**2+b_**2-blk)==0
    q=sp.cancel(br(F[j],F[k])/blk); assert q.is_Rational and q!=0; cof[j,k]=q
    for m in (j,k):
        if sp.rem(nrm(F[m]),blk,T)==0: iso2.add(((j,k),m))
    if not ((j,k)==(1,2) or j==0):
        for m in (j,k): assert sp.rem(nrm(F[m]),blk,T)!=0
print('T-chart isotropic incidences:',sorted(iso2))
print('cofactors',cof)
# norms vs n_j : degree of defect
for j in range(5):
    nj=sp.expand(sp.prod([((T+1)**2+h[e]**2*T**2)/(1+h[e]**2) for e in E if j in e]))
    dl=sp.Poly(sp.expand(nrm(F[j])-nj),T)
    print('mark',j+1,'deg |F|^2 =',sp.Poly(nrm(F[j]),T).degree(),'defect degree',dl.degree() if not dl.is_zero else 'exact')
import pickle; pickle.dump((F,h,E,cof),open('five_F.pkl','wb'))
