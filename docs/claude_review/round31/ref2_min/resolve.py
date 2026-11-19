# Resultant-based analysis for the classes where Groebner saturation was too slow.
import sys, itertools, sympy as sp, mpmath as mp
from probe import build, out, a, b
mp.mp.dps=60
def analyse(n):
    H,Y,V,rows=build(n)
    M=sp.Matrix([r[:3] for r in rows]); cc=sp.Matrix([r[3] for r in rows])
    s,t=sp.symbols('s t')
    Dpairs=list(itertools.combinations(range(8),2))
    print('class',n,'nq',out[n]['nq'],'nmid',out[n]['nm'],'H',[tuple(H[i,:]) for i in range(8)])
    for R in itertools.combinations(range(5),3):
        MR=M.extract(list(R),[0,1,2]); dR=sp.expand(MR.det())
        if dR==0: continue
        Xn=[sp.expand(x) for x in -(MR.adjugate()*cc.extract(list(R),[0]))]
        conds=[sp.expand(sum(M[k,q]*Xn[q] for q in range(3))+cc[k]*dR) for k in range(5) if k not in R]+[sp.expand(Xn[0]*Xn[2]-Xn[1]**2)]
        badf=[f for f,_ in sp.factor_list(sp.expand(dR*sp.Mul(*Y)))[1] if f.free_symbols]
        changed=True
        while changed:
            changed=False; g=sp.gcd_list(conds)
            for f in badf:
                q,r=sp.div(g,f,a,b)
                if r==0: conds=[sp.quo(cnd,f,a,b) for cnd in conds]; changed=True; break
        g=sp.gcd_list(conds)
        X=[x/dR for x in Xn]
        print(' R',R,'dR=',sp.factor(dR))
        # curve components
        if sp.Poly(g,a,b).total_degree()>0:
            for f,_ in sp.factor_list(g)[1]:
                if not f.free_symbols: continue
                print('   curve component',f)
                sol=sp.solve(f,a) if f.has(a) else None
                if sol and len(sol)==1:
                    sub={a:sol[0]}
                    Ys=[sp.factor(y.subs(sub)) for y in Y]; Xs=[sp.factor(x.subs(sub)) for x in X]
                    print('     on curve: Y=',Ys); print('     X=',Xs)
                    zeroD=[p for p in Dpairs if sp.simplify(((Y[p[1]]*V[p[0]]-Y[p[0]]*V[p[1]])**2+(Y[p[0]]-Y[p[1]])**2*Y[p[0]]**2*Y[p[1]]**2).subs({s**2:0}).subs(sub))==0]
                    # exact: D_ij as linear form in X along curve
                    Dl=[]
                    for (i,j) in Dpairs:
                        e=sp.Poly(sp.expand((Y[j]*V[i]-Y[i]*V[j])**2),s,t); val=0
                        for (p_,q_),c in e.terms(): val+=c*{(2,0):X[0],(1,1):X[1],(0,2):X[2]}[(p_,q_)]
                        val+= (Y[i]-Y[j])**2*Y[i]**2*Y[j]**2
                        Dl.append(sp.factor(sp.cancel(val.subs(sub))))
                    zero=[Dpairs[k] for k,d in enumerate(Dl) if d==0]
                    print('     pairs with D identically 0 on the curve:',zero)
                else: print('     (nonlinear component, manual)')
        cq=[sp.quo(cnd,g,a,b) for cnd in conds]
        cq=[c for c in cq if c!=0]
        r1=sp.resultant(cq[0],cq[1],a); r2=sp.resultant(cq[0],cq[2],a); r=sp.gcd(r1,r2)
        if r==0: print('   resultant identically zero?'); continue
        roots=sp.Poly(r,b).real_roots() if sp.Poly(r,b).degree()>0 else []
        print('   finite part: b-candidates (real):',len(set(roots)))
        for b0 in sorted(set(roots),key=lambda x: float(x)):
            bn=mp.mpf(str(sp.N(b0,70)))
            pa=sp.Poly(cq[0].subs(b,sp.Float(str(bn),70)),a)
            if pa.degree()<=0: print('    b0',sp.N(b0,8),'cq0 independent of a; manual'); continue
            for an in mp.polyroots([mp.mpf(str(c)) if c.is_real else mp.mpc(str(sp.re(c)),str(sp.im(c))) for c in pa.all_coeffs()],maxsteps=400,extraprec=400):
                if abs(mp.im(an))>1e-25: continue
                an=mp.re(an); sub={a:an,b:bn}
                res=max(abs(mp.mpf(str(sp.N(c.subs({a:sp.Float(str(an),70),b:sp.Float(str(bn),70)}),60)))) for c in cq)
                if res>1e-20: continue
                Yv=[mp.mpf(str(sp.N(y.subs({a:sp.Float(str(an),70),b:sp.Float(str(bn),70)}),60))) for y in Y]
                dv=sp.N(dR.subs({a:sp.Float(str(an),70),b:sp.Float(str(bn),70)}),60)
                Xv=[sp.N(x.subs({a:sp.Float(str(an),70),b:sp.Float(str(bn),70)}),40) if abs(dv)>1e-30 else None for x in X]
                print('    real common zero a=%s b=%s minY=%s X=%s'%(mp.nstr(an,12),mp.nstr(bn,12),mp.nstr(min(Yv),8),[sp.N(x,8) if x is not None else None for x in Xv]))
if __name__=='__main__':
    for n in map(int,sys.argv[1].split(',')): analyse(n); sys.stdout.flush()
