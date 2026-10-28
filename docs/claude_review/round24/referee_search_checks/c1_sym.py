# Independent symbolic check of Theorem A steps (a)-(c) (sufficiency direction).
import sympy as sp, pickle, time
b0,b1,b2,t,X,I = sp.symbols('b0 b1 b2 t X I')
b3 = -(b0+b1+b2)
bs=[b0,b1,b2,b3]
T2 = b0*b1*b2*b3
def red(e):
    # reduce I^2=-1 and t^2=T2
    e=sp.expand(e)
    p=sp.Poly(e, I, t)
    out=0
    for (ei,et),c in p.terms():
        out += c * (sp.Integer(-1))**(ei//2) * I**(ei%2) * T2**(et//2) * t**(et%2)
    return sp.expand(out)
a=[t/b for b in bs]
p=[a[r]+I*bs[r] for r in range(4)]
pc=[a[r]-I*bs[r] for r in range(4)]
def ReIm(e):
    e=red(e); e=sp.expand(e)
    return e.coeff(I,0), e.coeff(I,1)
lab={'123':p[0],'4':-pc[0],'1':p[1],'2':p[2],'3':p[3],'23':-pc[1],'13':-pc[2],'12':-pc[3]}
for j in '123':
    P=1
    for S,r in lab.items():
        if j in S: P*= (X-r)
    re,im=ReIm(P)
    im=sp.together(im); num,den=sp.fraction(im)
    pn=sp.Poly(sp.expand(num),X)
    print('Im P_%s : X-degree of numerator'%j, pn.degree(), ' const term nonzero:', pn.coeff_monomial(1)!=0)
n=1
for r in range(4): n*= (X-p[r])*(X+pc[r])
ren,imn=ReIm(n)
imn_p=sp.Poly(sp.together(imn).as_numer_denom()[0],X)
print('Im n numerator monomials (X-degrees):',[m[0] for m in imn_p.monoms()])
ren_p=sp.Poly(sp.together(ren).as_numer_denom()[0],X)
print('Re n numerator X-degrees:',[m[0] for m in ren_p.monoms()])
sig=-sp.prod([red(p[r]*pc[r]) for r in range(4)])
sig=red(sig)
print('sigma has t?', sig.has(t))
L=sp.together(ren - imn**2/(4*sig) + sig)
Ln,Ld=sp.fraction(L)
Ln=red(Ln); Ld=red(Ld)
print('L numerator has t?',Ln.has(t),' denominator has t?',Ld.has(t))
LP=sp.Poly(Ln,X)
print('L numerator X-degrees:',sorted(set(m[0] for m in LP.monoms())))
pickle.dump((Ln,Ld),open('c1_L.pkl','wb'))
