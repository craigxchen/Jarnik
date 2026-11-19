"""(c) two Walsh cores on labels x and x^3 in F_{2^d}: rank, word weights, square count, e>=4d-2."""
import sys
sys.path.insert(0,'.')
from gf2 import mul, mod, deg, deriv, is_irreducible
def field(d):
    # find irreducible modulus of degree d (own test)
    for f in range(1<<d, 1<<(d+1)):
        if f & 1 and is_irreducible(f): return f
for d in (3,5,7):
    f=field(d); m=1<<d
    fm=lambda a,b: mod(mul(a,b),f)
    cube=[fm(fm(x,x),x) for x in range(m)]
    assert len(set(cube))==m
    tr=lambda a: (sum(1 for _ in [0]) and None)
    def Tr(a):
        s=0; x=a
        for _ in range(d):
            s^=x; x=fm(x,x)
        assert s in (0,1); return s
    cols=[(0,u) for u in range(1,m)]+[(1,u) for u in range(1,m)]
    rows=[]
    for x in range(m):
        row=0
        for k,(c,u) in enumerate(cols):
            b=Tr(fm(u, x if c==0 else cube[x]))
            row|=b<<k
        rows.append(row)
    piv={}
    for r in rows[1:]:
        r^=rows[0]
        while r:
            l=r.bit_length()
            if l not in piv: piv[l]=r;break
            r^=piv[l]
    rk=len(piv)
    # enumerate codewords of the anchored code: words = combos of rows
    basis=list(piv.values())
    ws=set()
    for g in range(1,1<<rk):
        w=0
        for i in range(rk):
            if g>>i&1: w^=basis[i]
        a=bin(w & ((1<<(m-1))-1)).count('1'); b=bin(w>>(m-1)).count('1')
        ws.add((a,b))
    print(f"d={d} m={m}: anchored rank={rk} (2d={2*d}); (weight on core1, weight on core2) over nonzero words: {sorted(ws)}")
    # polarisation identity check: Tr(B x^2 y)+Tr(B x y^2) = Tr(x*((By)^(2^(d-1)) + B y^2))
    import random
    R=random.Random(d)
    for _ in range(200):
        B,x,y=R.randrange(m),R.randrange(m),R.randrange(m)
        lhs=Tr(fm(B,fm(fm(x,x),y)))^Tr(fm(B,fm(x,fm(y,y))))
        z=fm(B,y)
        for _ in range(d-1): z=fm(z,z)
        rhs=Tr(fm(x, z ^ fm(B,fm(y,y))))
        assert lhs==rhs
print("polarisation identity ok")
# number of nonzero squares of degree <= e/2-1 in F_2[t] vs 2^ceil(e/4)-1
for e in range(2,23,2):
    cnt=sum(1 for g in range(1,1<<12) if deg(mul(g,g))<=e//2-1)
    assert cnt==2**(-(-e//4))-1, (e,cnt)
print("square counts ok: #nonzero squares of degree <= e/2-1 = 2^ceil(e/4)-1 for even e<=22")
for d in (3,5,7,9):
    m=2**d
    e=next(e for e in range(2,200,2) if 2**(-(-e//4))-1>=m-1)
    print(f"d={d}: least even e with 2^ceil(e/4)-1>=m-1 is {e} (4d-2={4*d-2}); W>=2(m-1)e={2*(m-1)*e}, m log2 m={m*d}")
