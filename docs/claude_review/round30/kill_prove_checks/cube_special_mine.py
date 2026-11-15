# Independent elimination: real y on {sum y=0, prod y=1} with an extra coincidence d (not in Z(1,1,1,1)),
# sum d y=0, sum d g=0, g=3/y-y^3. Uses resultants (not Groebner).
import itertools, sympy as sp, mpmath as mp, flint
flint.ctx.prec=200
mp.mp.dps=60
u,v=sp.symbols('u v')
seen=set(); reps=[]
for d in itertools.product(range(-2,3),repeat=4):
    if len(set(d))==1: continue
    # normalise: shift so min entry 0? conditions invariant under d->d+lam*1 and d->-d and S_4
    best=None
    for lam in range(-4,5):
        for sg in (1,-1):
            e=tuple(sorted(sg*(x+lam) for x in d))
            if best is None or e<best: best=e
    if best in seen: continue
    seen.add(best); reps.append(d)
print('d classes',len(reps))
special=[]
for d in reps:
    # parametrise the plane sum y=0, sum d y=0 by (u,v)
    Y=sp.symbols('y0:4')
    i,j=next((i,j) for i in range(4) for j in range(i+1,4) if d[i]!=d[j])
    fr=[k for k in range(4) if k not in (i,j)]
    sol=sp.solve([sum(Y), sum(di*yi for di,yi in zip(d,Y))],[Y[i],Y[j]],dict=True)[0]
    yy=[sp.expand(sp.sympify(Yk).subs(sol).subs({Y[fr[0]]:u,Y[fr[1]]:v})) for Yk in Y]
    f1=sp.expand(yy[0]*yy[1]*yy[2]*yy[3]-1)
    g=sum(di*(3/yk-yk**3) for di,yk in zip(d,yy))
    f2=sp.numer(sp.together(g)); f2=sp.expand(f2)
    if f2==0: print('identically zero for',d); continue
    G=sp.gcd(f1,f2)
    if sp.Poly(G,u,v).total_degree()>0:
        # curve component: sample it and test degeneracy
        degs=[]
        for uu in [sp.Rational(1,3),sp.Rational(7,5),sp.Rational(-2,3),sp.Rational(5,2)]:
            for vv in sp.Poly(G.subs(u,uu),v).nroots(n=30):
                yv=[complex(sp.N(yk.subs({u:uu,v:vv}))) for yk in yy]
                degs.append(min(abs(yv[a]**2-yv[b]**2) for a in range(4) for b in range(a+1,4)))
        print('curve component for d=',d,'gcd=',sp.factor(G),' max over samples of min|y_a^2-y_b^2|:',max(degs) if degs else None)
        f1=sp.cancel(f1/G); f2=sp.cancel(f2/G)
    R=sp.Poly(sp.resultant(f1,f2,v),u)
    if R.is_zero: print('resultant still zero',d); continue
    for fac,_ in sp.factor_list(R.as_expr())[1]:
        P=sp.Poly(fac,u)
        if P.degree()<1: continue
        fp=flint.fmpz_poly([int(c) for c in reversed(P.all_coeffs())])
        for r,_m in fp.complex_roots():
            if not r.imag.contains(0): continue   # rigorous: nonreal roots excluded
            ur=sp.Float(r.real.mid().str(50,radius=False),50)
            p1=sp.Poly(sp.expand(f1.subs(u,ur)),v)
            for vr in sp.Poly(p1,v).nroots(n=45,maxsteps=2000):
                if abs(sp.im(vr))>1e-25: continue
                vr=sp.re(vr)
                if abs(f2.subs({u:ur,v:vr}))>1e-15: continue
                yv=[mp.mpf(str(sp.N(yk.subs({u:ur,v:vr}),45))) for yk in yy]
                if min(abs(x) for x in yv)<1e-20: continue
                ay=sorted(abs(x) for x in yv)
                nondeg=min(ay[i+1]-ay[i] for i in range(3))>1e-20
                special.append((d,yv,nondeg))
print('real solutions found',len(special),' nondegenerate:',sum(1 for s in special if s[2]))
sets=set()
for d,yv,nd in special:
    if nd: sets.add(tuple(mp.nstr(x,12) for x in sorted(abs(t) for t in yv)))
print('distinct |y| sets (nondeg):',sets)
import pickle; pickle.dump([(d,[str(x) for x in yv],nd) for d,yv,nd in special],open('cube_special_mine.pkl','wb'))
