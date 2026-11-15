# All base points of a cube N: differences d=alpha-alpha' in {-2..2}^4. Check no nondegenerate solution
# of sum d y=0, sum d g=0 (g=3/y-y^3) on the surface unless d in Z(1,1,1,1).
import itertools, sympy as sp
y0,y1,y2,u=sp.symbols('y0 y1 y2 u')
y3=-(y0+y1+y2); Y=[y0,y1,y2,y3]
surf=sp.expand(y0*y1*y2*y3-1)
nd=1
for i in range(4):
    for j in range(i+1,4): nd*=(Y[i]**2-Y[j]**2)
sat=sp.expand(u*nd-1)
seen=set(); bad=[]
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
    # univariate in y2?
    uni=[g for g in el if g.free_symbols<= {y2}]
    realroots=None
    if uni:
        realroots=sp.Poly(uni[0],y2).count_roots()  # real roots count
    bad.append((key,realroots))
    print(key,'solutions exist; univariate y2 real roots:',realroots, sp.factor(uni[0]) if uni else el[-1],flush=True)
print('done; classes with complex nondeg solutions:',len(bad))
