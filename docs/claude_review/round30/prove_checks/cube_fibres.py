# For generic and for every real special member (from d in {-2..2}^4 coincidences), compute the max fibre
# of s -> moments over all 256 points of N (high precision). Fibre = cluster with any base point.
import itertools, sympy as sp, mpmath as mp, random
mp.mp.dps=60
y0,y1,y2,u=sp.symbols('y0 y1 y2 u')
y3=-(y0+y1+y2); Y=[y0,y1,y2,y3]
surf=sp.expand(y0*y1*y2*y3-1)
nd=1
for i in range(4):
    for j in range(i+1,4): nd*=(Y[i]**2-Y[j]**2)
sat=sp.expand(u*nd-1)
eps=[1,-1,-1,1]
def maxfibre(yv):
    # roots: (k,+): r_k, (k,-): -r_k ; r_k=eps/y+iy. point picks root or conj for each of 8 positions
    roots=[]
    for k in range(4):
        r=mp.mpc(eps[k]/yv[k],yv[k]); roots+= [r,-r]
    from collections import Counter
    C=Counter()
    for s in itertools.product([1,-1],repeat=8):
        v=tuple(mp.nstr(sum(s[j]*mp.im(roots[j]**m) for j in range(8)),25) for m in (1,2,3))
        v=tuple('0.0' if x in('-0.0',) else x for x in v)
        C[v]+=1
    return max(C.values()), sorted(C.values())[-3:]
# generic
random.seed(1)
for trial in range(3):
    # pick y0,y1 random, solve y2 from product=1 (cubic in y2)
    while True:
        a=mp.mpf(random.uniform(-3,3)); b=mp.mpf(random.uniform(-3,3))
        # y2*(-(a+b+y2))*a*b=1 -> a b y2^2 + a b (a+b) y2 + 1 = 0
        disc=(a*b*(a+b))**2-4*a*b
        if disc>0: break
    c=(-a*b*(a+b)+mp.sqrt(disc))/(2*a*b); yv=[a,b,c,-(a+b+c)]
    print('generic',[mp.nstr(x,6) for x in yv],maxfibre(yv),flush=True)
seen=set(); pts=[]
for d in itertools.product(range(-2,3),repeat=4):
    if all(x==0 for x in d) or len(set(d))==1: continue
    key=min(d,tuple(-x for x in d))
    if key in seen: continue
    seen.add(key)
    e1=sum(d[k]*Y[k] for k in range(4))
    n3=sp.numer(sp.together(sum(d[k]*(3/Y[k]-Y[k]**3) for k in range(4))))
    sols=sp.solve_poly_system([sp.expand(e1),sp.expand(n3),surf],y0,y1,y2) if False else None
    G=sp.groebner([sp.expand(e1),sp.expand(n3),surf,sat],u,y0,y1,y2,order='lex')
    if G.exprs==[1]: continue
    el=[g for g in G.exprs if not g.has(u)]
    uni=[g for g in el if g.free_symbols<={y2}][0]
    for rt in sp.Poly(uni,y2).real_roots():
        v2=sp.N(rt,70)
        # back-substitute: solve remaining eqs numerically in y0,y1
        rest=[sp.expand(e1),sp.expand(n3),surf]
        eqs=[e.subs(y2,v2) for e in rest]
        # solve e1 for y0 (linear if d0!=0) else use groebner elements
        cand=[]
        G2=sp.groebner([e.subs(y2,rt) for e in rest]+[],y0,y1,order='lex')
        try:
            s=sp.solve([g for g in el if g.has(y0) or g.has(y1)] and [g.subs(y2,rt) for g in el if g.has(y0) or g.has(y1)],[y0,y1],dict=True)
        except Exception as ex:
            s=[]
        for so in s:
            vals=[sp.N(so[y0],60),sp.N(so[y1],60),v2]
            if any(abs(sp.im(x))>1e-30 for x in vals): continue
            yv=[mp.mpf(str(sp.re(x))) for x in vals]; yv.append(-(yv[0]+yv[1]+yv[2]))
            if abs(yv[0]*yv[1]*yv[2]*yv[3]-1)>1e-40: continue
            mf=maxfibre(yv)
            pts.append((key,yv,mf))
            print(key,[mp.nstr(x,8) for x in yv],'maxfibre',mf,flush=True)
print('special real points:',len(pts),'max fibre overall',max([p[2][0] for p in pts]+[0]))
