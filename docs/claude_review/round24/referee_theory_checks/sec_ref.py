# Independent check of the building blocks in Prop 3.1(b): support in S_m, iota-parity, order at e,
# and of the full table for m = 1..9 by exact products (orders computed directly, not by additivity).
import itertools, sympy as sp
from collections import defaultdict
def mul(f,g):
    h=defaultdict(int)
    for s,a in f.items():
        for t,b in g.items():
            h[tuple(x+y for x,y in zip(s,t))]+=a*b
    return {k:v for k,v in h.items() if v}
def e(j): return tuple(1 if i==j else 0 for i in range(5))
def lin(d): return {k:v for k,v in d.items() if v}
def diffmono(a,b): # x_a - x_b as Laurent poly (not in M; combine later)
    return {e(a):1, e(b):-1}
def shift(f,s): return {tuple(x+y for x,y in zip(k,s)):v for k,v in f.items()}
def Vabc(a,b,c):
    f=mul(mul(diffmono(a,b),diffmono(b,c)),diffmono(a,c))
    return shift(f,tuple(-1 if i in (a,b,c) else 0 for i in range(5)))
def D(a,b):
    s=tuple((1 if i==a else -1 if i==b else 0) for i in range(5))
    return {s:1, tuple(-x for x in s):1, (0,)*5:-2}
V={(0,)*5:1}
for a,b in itertools.combinations(range(5),2): V=mul(V,diffmono(a,b))
V=shift(V,(-2,)*5)
y=sp.symbols('y1:5')
def order(f,maxd=40):
    # ord_e f = min j with sum c_s <s,y>^j != 0, y5=0
    for j in range(maxd):
        tot=0
        for s,c in f.items():
            tot+=c*sum(si*yi for si,yi in zip(s[:4],y))**j
        if sp.expand(tot)!=0: return j
    return None
def inS(f,m): return all(sum(s)==0 and sum(abs(x) for x in s)<=2*m for s in f)
def parity(f):
    g={tuple(-x for x in s):c for s,c in f.items()}
    if g==f: return +1
    if {k:-v for k,v in g.items()}==f: return -1
    return 0
def mmin(f):  # least m with f in S_m
    return max((sum(abs(x) for x in s)+1)//2 for s in f)
blocks={'V':V,'V123':Vabc(0,1,2),'V145':Vabc(0,3,4),'V245':Vabc(1,3,4),'D12':D(0,1),'D45':D(3,4)}
for n,f in blocks.items():
    print(n,'min m',mmin(f),'parity',parity(f),'ord',order(f),flush=True)
def power(f,k):
    r={(0,)*5:1}
    for _ in range(k): r=mul(r,f)
    return r
B=blocks
for m in range(1,8):
    s,r=divmod(m,3)
    if r==0: ev=power(V,s); od=mul(power(V,s-1),mul(B['V123'],mul(B['V145'],B['V245'])))
    elif r==1: ev=mul(power(V,s),B['D12']); od=mul(power(V,s),B['V123'])
    else: ev=mul(power(V,s),mul(B['V123'],B['V145'])); od=mul(power(V,s),mul(B['V123'],B['D45']))
    oe=order(ev,60) if m<=4 else None; oo=order(od,60) if m<=4 else None
    print(f"m={m} t={(10*m)//3} even: inS={inS(ev,m)} par={parity(ev)} ord={oe} | odd: inS={inS(od,m)} par={parity(od)} ord={oo}",flush=True)
