# Special members of the cube locus: extra zero-sum sets.
# zero-sum F <-> (a,b): a_k=c+ - c-, b_k = c+ + c- ; conditions sum b eps=0, sum a y=0, sum a(3/y-y^3)=0
import itertools, sympy as sp
y0,y1,y2=sp.symbols('y0 y1 y2')
y3=-(y0+y1+y2)
Y=[y0,y1,y2,y3]
eps=[1,-1,-1,1]
surf=sp.expand(y0*y1*y2*y3-1)
seen=set(); results=[]
for a in itertools.product([-1,0,1],repeat=4):
    if all(x==0 for x in a): continue
    if len(set(a))==1: continue   # +-(1,1,1,1)
    key=min(a,tuple(-x for x in a))
    if key in seen: continue
    seen.add(key)
    # does some b exist? b_k=1 where a_k!=0, b_k in{0,2} where a_k=0
    zeros=[k for k in range(4) if a[k]==0]
    okb=False
    for bz in itertools.product([0,2],repeat=len(zeros)):
        b=[1 if a[k]!=0 else None for k in range(4)]
        for k,v in zip(zeros,bz): b[k]=v
        if sum(b[k]*eps[k] for k in range(4))==0: okb=True
    e1=sum(a[k]*Y[k] for k in range(4))
    e3=sp.together(sum(a[k]*(3/Y[k]-Y[k]**3) for k in range(4)))
    n3=sp.numer(e3)
    G=sp.groebner([sp.expand(e1),sp.expand(n3),surf],y0,y1,y2,order='lex')
    results.append((a,okb,G.exprs[-1] if G.exprs!=[1] else 1))
    print(a,'b-ok' if okb else 'no-b', 'GB last:',sp.factor(G.exprs[-1]) if G.exprs!=[1] else 'INCONSISTENT', flush=True)
