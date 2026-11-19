# Judge's independent checks for Part V (Question B thread). sympy only.
from sympy import *
T=symbols('T',real=True); x=symbols('x')
def br(a,b):  # det of (Re a, Im a),(Re b, Im b) for polys with Gaussian-rational coeffs, T real
    a=expand(a); b=expand(b)
    ar,ai=a.as_real_imag(); br_,bi=b.as_real_imag()
    return expand(ar*bi-ai*br_)
def cnj(p): return expand(p.subs(I,-I))
# --- four-star, parametrised shape
v,t,k=Rational(3,2),2,2
f=(1-t**2)*(t**4-3*t**2+1); X=(k+f/k)/2; Y=(f/k-k)/(2*t)
w=[1]+[v/2*(1+X),v/2*(1-X),v/2*(-t**2+Y),v/2*(-t**2-Y)]
w=[nsimplify(z) for z in w]
s=[sum(prod(c) for c in __import__('itertools').combinations(w[1:],r)) for r in range(5)]
print('w',w,' s3-s1s2 =',s[3]-s[1]*s[2],' s1s4',s[1]*s[4])
G=expand(prod(T-I*z for z in w[1:])); A=-T/s[1]+I
u=[z**2 for z in w]
F=[expand(G+(ui-1)/ui*(ui+T**2)*A) for ui in u]
print('[G,A]=',simplify(br(G,A)))
ok=all(simplify(br(F[i],F[j])-s[4]*(1/u[i]-1/u[j])*(T**2+u[i]*u[j]))==0 for i in range(5) for j in range(i+1,5))
print('ten minor identities:',ok)
J={(i,j):T**2+u[i]*u[j] for i in range(5) for j in range(i+1,5)}
n=[expand(prod(J[e] for e in J if i in e)) for i in range(5)]
nrm=[expand(Fi*cnj(Fi)) for Fi in F]
for i in range(5):
    d=expand(nrm[i]-n[i]); print('mark',i+1,'deg n',degree(n[i],T),'deg defect',degree(d,T) if d!=0 else 'exact')
iso=[(e,i) for e in J for i in e if rem(nrm[i],J[e],T)==0]
print('isotropic incidences:',len(iso),sorted(set(e for e,_ in iso)))
# --- five-block curve
z=[Rational(3,2),3,Rational(1,3),Rational(-32,15),Rational(-5,54)]
u=[q**2 for q in z]; h=lambda j,kk:z[j]*z[kk]
G=expand(prod(x+I*h(0,j) for j in range(1,5)))
B=Rational(4500480,1226641)*x**2+Rational(1413120,1226641)+I*(x**2+Rational(26758,283071)*x+1)
def brx(a,b):
    a=expand(a);b=expand(b)
    ar,ai=a.as_real_imag(); b1,b2=b.as_real_imag(); return expand(ar*b2-ai*b1)
x=symbols('x',real=True); G=G.subs(symbols('x'),x); B=B.subs(symbols('x'),x)
print('[G,B]-(x-1)^6 =',simplify(brx(G,B)-(x-1)**6))
d=-G.subs(x,I)/(B.subs(x,I)*(u[1]-u[0])*(u[0]-u[2])); d=nsimplify(simplify(d))
print('d =',d,' equals claim:',d==Rational(-418284581,1009043070))
L=[expand(G+d*B*(uj-u[0])*(u[0]+x**2/uj)) for uj in u]
print('L2(i),L3(i):',simplify(L[1].subs(x,I)),simplify(L[2].subs(x,I)))
F=[expand(T**4*Lj.subs(x,1+1/T)/Lj.subs(x,1)) for Lj in L]
F=[expand(cancel(Fi)) for Fi in F]
Jt={(a,b):expand(((T+1)**2+h(a,b)**2*T**2)/(1+h(a,b)**2)) for a in range(5) for b in range(a+1,5)}
cof={}
okm=True
for (a,b),Je in Jt.items():
    q,r=div(br(F[a],F[b]),Je,T); okm&= (r==0 and degree(q,T)==0 and q!=0); cof[(a+1,b+1)]=q
print('five-block: ten minors = nonzero const * J:',okm)
nrm=[expand(Fi*cnj(Fi)) for Fi in F]
iso=sorted(set((a+1,b+1) for (a,b),Je in Jt.items() for i in (a,b) if rem(nrm[i],Je,T)==0))
inc=sum(1 for (a,b),Je in Jt.items() for i in (a,b) if rem(nrm[i],Je,T)==0)
print('isotropic blocks',iso,'incidences',inc)
n=[expand(prod(Jt[e] for e in Jt if i in e)) for i in range(5)]
for i in range(5):
    dd=expand(nrm[i]-n[i]); print('mark',i+1,'deg defect',degree(dd,T) if dd!=0 else 'exact', 'lead ratio', LC(nrm[i],T)/LC(n[i],T))
