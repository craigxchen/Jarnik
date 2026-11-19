# Judge's reproduction of the bounds referee's 8-point cluster with nonconstant translation B:
# V = F_2[t]_{<4}, B = t^2+t, Y = t^11 + theta*B, all 8 forms l with l(B)=1, every choice of sides.
import itertools
MUL=[[0,0,0,0],[0,1,2,3],[0,2,3,1],[0,3,1,2]]
INV=[0,1,3,2]
CONJ=[0,1,3,2]
def trim(a):
    while a and a[-1]==0: a.pop()
    return a
def pmul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y: r[i+j]^=MUL[x][y]
    return trim(r)
def padd(a,b):
    r=[0]*max(len(a),len(b))
    for i,x in enumerate(a): r[i]^=x
    for i,x in enumerate(b): r[i]^=x
    return trim(r)
def pmod(a,b):
    a=a[:]; db=len(b)-1; il=INV[b[-1]]
    while len(a)-1>=db and a:
        c=MUL[a[-1]][il]; s=len(a)-1-db
        for i,y in enumerate(b):
            if y: a[s+i]^=MUL[c][y]
        trim(a)
    return a
def pgcd(a,b):
    while b: a,b=b,pmod(a,b)
    il=INV[a[-1]]; return [MUL[il][x] for x in a]
def pdiv(a,b):
    a=a[:]; db=len(b)-1; q=[0]*(len(a)-db); il=INV[b[-1]]
    while len(a)-1>=db and a:
        c=MUL[a[-1]][il]; s=len(a)-1-db; q[s]=c
        for i,y in enumerate(b):
            if y: a[s+i]^=MUL[c][y]
        trim(a)
    assert not a; return q
def conj(a): return [CONJ[x] for x in a]
Lexp=11; B=0b110
V=list(range(16))          # bit i = coefficient of t^i
def bits(v,n): return [(v>>i)&1 for i in range(n)]
def Yv(v):                 # t^L + theta*B + v
    c=[0]*(Lexp+1); c[Lexp]=1
    for i in range(4):
        if (v>>i)&1: c[i]^=1
        if (B>>i)&1: c[i]^=2
    return c
forms=[l for l in range(1,16) if bin(l&B).count('1')%2==1]
assert len(forms)==8
def side(l,s): return [v for v in V if bin(l&v).count('1')%2==s]
F={}
for l in forms:
    for s in (0,1):
        p=[1]
        for v in side(l,s): p=pmul(p,Yv(v))
        F[(l,s)]=p
PV=[1]
for v in V: PV=pmul(PV,Yv(v))
assert all(pmul(F[(l,0)],conj(F[(l,0)]))==PV for l in forms)   # conj swaps sides: common norm P_V(Y)
best=None
for sides in itertools.product((0,1),repeat=8):
    if sides[0]==1: continue      # global conjugation symmetry
    pts=[F[(l,s)] for l,s in zip(forms,sides)]
    G=pts[0]
    for p in pts[1:]: G=pgcd(G,p)
    Z=[pdiv(p,G) for p in pts]
    D=len(Z[0])-1
    md=max(len(padd(Z[i],Z[j]))-1 for i in range(8) for j in range(i+1,8))
    if md<=D/2 and (best is None or D<best[0]): best=(D,sides,len(G)-1,md)
print("best 8-point c=0 cluster: D=%s sides=%s gcd deg=%s maxdiff=%s"%best)
D,sides,_,_=best
pts=[pdiv(F[(l,s)],G) for l,s in zip(forms,sides)] if False else None
pts=[F[(l,s)] for l,s in zip(forms,sides)]
G=pts[0]
for p in pts[1:]: G=pgcd(G,p)
Z=[pdiv(p,G) for p in pts]
Ns=[pmul(z,conj(z)) for z in Z]
print("reduced: distinct",len(set(map(tuple,Z)))==8,"monic",all(z[-1]==1 for z in Z),
      "norms equal",all(n==Ns[0] for n in Ns),"norm in F_2[t]",all(x in (0,1) for x in Ns[0]),"deg N",len(Ns[0])-1)
print("pair diffs",sorted(len(padd(Z[i],Z[j]))-1 for i in range(8) for j in range(i+1,8)))
