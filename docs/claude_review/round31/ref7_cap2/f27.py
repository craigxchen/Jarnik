# own F_27 = F_3[a]/(a^3+2a+1); elements as ints 0..26 (base-3 digits c0+3c1+9c2)
import itertools
def to(v): return (v%3, (v//3)%3, (v//9)%3)
def fr(c): return c[0]%3+3*(c[1]%3)+9*(c[2]%3)
def add(x,y): a,b=to(x),to(y); return fr([a[i]+b[i] for i in range(3)])
def neg(x): a=to(x); return fr([-a[i] for i in range(3)])
def mul(x,y):
    a,b=to(x),to(y); c=[0]*5
    for i in range(3):
        for j in range(3): c[i+j]+=a[i]*b[j]
    # a^3 = -2a-1 = a+2 ; a^4 = a^2+2a
    for d in (4,3):
        t=c[d]; c[d]=0; c[d-3]+=2*t; c[d-2]+=t   # a^d = a^{d-3}*(a+2)
    return fr(c[:3])
ADD=[[add(x,y) for y in range(27)] for x in range(27)]
MUL=[[mul(x,y) for y in range(27)] for x in range(27)]
def pw(x,e):
    r=1
    for _ in range(e): r=MUL[r][x]
    return r
# check field: multiplicative group cyclic of order 26
assert all(pw(x,26)==1 for x in range(1,27))
assert any(len({pw(x,e) for e in range(26)})==26 for x in range(1,27))
CHI=[0]+[1 if pw(x,13)==1 else -1 for x in range(1,27)]
assert CHI[neg(1)]==-1
