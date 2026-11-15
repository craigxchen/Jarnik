# Complete list of real nondegenerate special cube members (coincidences d in {-2..2}^4) and their max fibres.
import itertools, sympy as sp, mpmath as mp
from collections import Counter
mp.mp.dps=50
exec(open('cube_fibres2.py').read().split('random.seed(2)')[0])
y0,y1,y2,u=sp.symbols('y0 y1 y2 u')
y3=-(y0+y1+y2); Y=[y0,y1,y2,y3]
surf=sp.expand(y0*y1*y2*y3-1)
nd=1
for i in range(4):
    for j in range(i+1,4): nd*=(Y[i]**2-Y[j]**2)
sat=sp.expand(u*nd-1)
seen=set(); allpts={}
def eval_ok(polys,vals):
    return all(abs(sp.N(p.subs(vals),40))<1e-25 for p in polys)
for d in itertools.product(range(-2,3),repeat=4):
    if all(x==0 for x in d) or len(set(d))==1: continue
    key=min(d,tuple(-x for x in d))
    if key in seen: continue
    seen.add(key)
    e1=sum(d[k]*Y[k] for k in range(4))
    n3=sp.numer(sp.together(sum(d[k]*(3/Y[k]-Y[k]**3) for k in range(4))))
    G=sp.groebner([sp.expand(e1),sp.expand(n3),surf,sat],u,y0,y1,y2,order='lex')
    if G.exprs==[1]: continue
    el=[g for g in G.exprs if not g.has(u)]
    uni=[g for g in el if g.free_symbols<={y2}][0]
    for r2 in sp.Poly(uni,y2).nroots(n=40):
        if abs(sp.im(r2))>1e-20: continue
        r2=sp.re(r2)
        p1=[g.subs(y2,r2) for g in el if g.free_symbols<={y1,y2} and g.has(y1)]
        p1s=min(p1,key=lambda p:sp.degree(p,y1))
        for r1 in sp.Poly(sp.N(p1s,40),y1).nroots(n=40):
            if abs(sp.im(r1))>1e-20: continue
            r1=sp.re(r1)
            p0=[g.subs({y2:r2,y1:r1}) for g in el if g.has(y0)]
            p0s=min(p0,key=lambda p:sp.degree(p,y0))
            for r0 in sp.Poly(sp.N(p0s,40),y0).nroots(n=40):
                if abs(sp.im(r0))>1e-20: continue
                r0=sp.re(r0)
                vals={y0:r0,y1:r1,y2:r2}
                if not eval_ok([sp.expand(e1),sp.expand(n3),surf],vals): continue
                yv=[mp.mpf(str(r0)),mp.mpf(str(r1)),mp.mpf(str(r2))]; yv.append(-sum(yv))
                k=tuple(mp.nstr(v,12) for v in yv)
                if k not in allpts:
                    allpts[k]=fibres(yv)
                    print(key,k,allpts[k],flush=True)
print('TOTAL real special members',len(allpts),'abs-value sets',set(tuple(sorted(mp.nstr(abs(mp.mpf(x)),8) for x in k)) for k in allpts))
print('fibre profiles',Counter(tuple(v) for v in allpts.values()))
