# Independent check of Part V construction: F_2[t] polys as ints (bit i = coeff t^i);
# F_4[t] elements as pairs (A,B) meaning A + B*theta, theta^2=theta+1.
import itertools, sys
def mul2(a,b):
    r=0
    while b:
        if b&1: r^=a
        a<<=1; b>>=1
    return r
def deg(a): return a.bit_length()-1  # deg 0 = -1 (stand-in for -inf)
def mul4(x,y):
    A,B=x; C,D=y
    bd=mul2(B,D)
    return (mul2(A,C)^bd, mul2(A,D)^mul2(B,C)^bd)
def add4(x,y): return (x[0]^y[0],x[1]^y[1])
def conj(x): return (x[0]^x[1],x[1])
def deg4(x): return max(deg(x[0]),deg(x[1]))
def run(r):
    m=2**r; d=r+1; Q=2**d
    V=list(range(Q))   # polys of degree <d
    # functionals l: coefficient vector w in F_2^d, l(v)=popcount(v&w)%2, l(1)=w&1 must be 1
    forms=[w for w in range(Q) if w&1]
    assert len(forms)==m
    L=(d-2)*2**(d-1)
    Y=((1<<L),1)  # t^L + theta
    pts=[]
    for w in forms:
        H=[v for v in V if bin(v&w).count('1')%2==0]
        z=(1,0)
        for v in H:
            if v: z=mul4(z,add4(Y,(v,0)))
        pts.append(z)
    N=None; norms=set()
    for z in pts:
        n=mul4(z,conj(z)); assert n[1]==0; norms.add(n[0])
    assert len(norms)==1
    N=norms.pop(); D=deg4(pts[0])
    assert all(deg4(z)==D for z in pts) and deg(N)==2*D==(Q-2)*L
    mx=max(deg4(add4(a,b)) for a,b in itertools.combinations(pts,2))
    assert len(set(pts))==m
    # test minimality of L: try smaller L' 
    print(f"m={m} L={L} D={D} degN={deg(N)} maxdiff={mx} D/2={D/2} ok={mx<=D/2}")
    return D
if __name__=="__main__":
    for r in range(2,int(sys.argv[1])+1): run(r)
