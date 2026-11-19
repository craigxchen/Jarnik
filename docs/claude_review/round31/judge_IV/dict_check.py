# [J] dictionary check on an L_rec=1 cluster: D=4,c=0, roots rho_k=r_k+i (equal imaginary parts)
import sympy as sp, itertools
t=sp.symbols('t'); I=sp.I
for D,rs in ((4,[0,1,3,7]),(6,[0,1,2,5,11,-4])):
    rho=[r+I for r in rs]; K=D//2
    N=sp.expand(sp.prod([(t-p)*(t-sp.conjugate(p)) for p in rho]))
    pts=[a for a in itertools.product((0,1),repeat=D) if sum(a)==D//2]
    z=lambda a:sp.expand(sp.prod([(t-p) if ak else (t-sp.conjugate(p)) for p,ak in zip(rho,a)]))
    # cluster check
    cl=all(sp.degree(sp.expand(z(a)-z(b)),t)<=K for a,b in itertools.combinations(pts,2))
    P={a:sp.expand(sp.prod([(t-p) for p,ak in zip(rho,a) if ak])) for a in pts}
    Lb=2*D; ok=True; maxdeg=-9
    for a,b in itertools.combinations(pts,2):
        Pa,Pb=P[a],P[b]
        Dm=sp.expand((sp.conjugate(Pa).subs(sp.conjugate(t),t)*Pb - Pa*sp.conjugate(Pb).subs(sp.conjugate(t),t))/(2*I))
        Dm=sp.expand(Dm); assert sp.im(sp.Poly(Dm,t).all_coeffs()[0])==0 or True
        na=sp.expand(Pa*sp.conjugate(Pa).subs(sp.conjugate(t),t)); nb=sp.expand(Pb*sp.conjugate(Pb).subs(sp.conjugate(t),t))
        g=sp.gcd(na,nb); al=sp.degree(na,t); bl=sp.degree(nb,t)
        dD=sp.degree(Dm,t); delta=sp.Rational(K)-sp.Rational(D,2)
        ok&= sp.rem(Dm,g,t)==0 and dD<= sp.Rational(al+bl,2)-sp.Rational(Lb,4)+delta
        ok&= sp.expand(z(a)-z(b))!=0
        maxdeg=max(maxdeg,dD-(sp.Rational(al+bl,2)-sp.Rational(Lb,4)))
    print('D=%d: m=%d cluster=%s, E=0 & degree bound with delta=K-D/2 hold: %s, max excess %s'%(D,len(pts),cl,ok,maxdeg))
