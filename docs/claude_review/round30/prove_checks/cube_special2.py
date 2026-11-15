# Special members of the cube locus with |y_k| pairwise distinct (nondegenerate).
import itertools, sympy as sp
y0,y1,y2,u=sp.symbols('y0 y1 y2 u')
y3=-(y0+y1+y2)
Y=[y0,y1,y2,y3]
eps=[1,-1,-1,1]
surf=sp.expand(y0*y1*y2*y3-1)
nd=1
for i in range(4):
    for j in range(i+1,4): nd*= (Y[i]**2-Y[j]**2)
sat=sp.expand(u*nd-1)
seen=set()
for a in itertools.product([-1,0,1],repeat=4):
    if all(x==0 for x in a) or len(set(a))==1: continue
    key=min(a,tuple(-x for x in a))
    if key in seen: continue
    seen.add(key)
    zeros=[k for k in range(4) if a[k]==0]
    bs=[]
    for bz in itertools.product([0,2],repeat=len(zeros)):
        b=[1 if a[k]!=0 else None for k in range(4)]
        for k,v in zip(zeros,bz): b[k]=v
        if sum(b[k]*eps[k] for k in range(4))==0: bs.append(b)
    e1=sum(a[k]*Y[k] for k in range(4))
    n3=sp.numer(sp.together(sum(a[k]*(3/Y[k]-Y[k]**3) for k in range(4))))
    G=sp.groebner([sp.expand(e1),sp.expand(n3),surf,sat],u,y0,y1,y2,order='lex')
    if G.exprs==[1]:
        print(key,'nb=%d'%len(bs),'no nondegenerate solution',flush=True); continue
    el=[g for g in G.exprs if not g.has(u)]
    print(key,'nb=%d'%len(bs),'elim:',[sp.factor(g) for g in el][-2:],flush=True)
