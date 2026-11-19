# Independent check of the III.5 geometric pivot certificates: for the stated base point in F_{p^2},
# for all 36 matchings of the roots of V and Vbar, Res_b(E12,E23) != 0 (homogeneous resultants).
import itertools, sys
class Fq:
    p=None
    def __init__(s,x,y=0): s.x=x%Fq.p; s.y=y%Fq.p
    def __add__(s,o): o=o if isinstance(o,Fq) else Fq(o); return Fq(s.x+o.x,s.y+o.y)
    __radd__=__add__
    def __sub__(s,o): o=o if isinstance(o,Fq) else Fq(o); return Fq(s.x-o.x,s.y-o.y)
    def __rsub__(s,o): return Fq(o)-s
    def __neg__(s): return Fq(-s.x,-s.y)
    def __mul__(s,o): o=o if isinstance(o,Fq) else Fq(o); return Fq(s.x*o.x-s.y*o.y,s.x*o.y+s.y*o.x)
    __rmul__=__mul__
    def conj(s): return Fq(s.x,-s.y)
    def inv(s):
        n=(s.x*s.x+s.y*s.y)%Fq.p; ni=pow(n,-1,Fq.p); return Fq(s.x*ni,-s.y*ni)
    def __truediv__(s,o): return s*o.inv()
    def iszero(s): return s.x==0 and s.y==0
    def __pow__(s,k):
        r=Fq(1); b=s
        while k:
            if k&1: r=r*b
            b=b*b; k>>=1
        return r
    def __repr__(s): return f"{s.x}+{s.y}i"
def dec(n,p): return Fq(n%p,n//p)
def polmul(A,B):
    C=[Fq(0)]*(len(A)+len(B)-1)
    for i,a in enumerate(A):
        for j,b in enumerate(B): C[i+j]=C[i+j]+a*b
    return C
def lin(r): return [-r,Fq(1)]
def prod(rs):
    f=[Fq(1)]
    for r in rs: f=polmul(f,lin(r))
    return f
def det(M):
    n=len(M); M=[row[:] for row in M]; d=Fq(1)
    for c in range(n):
        piv=next((r for r in range(c,n) if not M[r][c].iszero()),None)
        if piv is None: return Fq(0)
        if piv!=c: M[c],M[piv]=M[piv],M[c]; d=-d
        d=d*M[c][c]; ic=M[c][c].inv()
        for r in range(c+1,n):
            f=M[r][c]*ic
            if not f.iszero():
                for k in range(c,n): M[r][k]=M[r][k]-f*M[c][k]
    return d
# polynomials in b as lists low->high
def padd(A,B):
    n=max(len(A),len(B)); return [ (A[i] if i<len(A) else Fq(0))+(B[i] if i<len(B) else Fq(0)) for i in range(n)]
def pneg(A): return [-a for a in A]
def W_coeffs(L,R):
    # N(b,X)=R(b)L(X)-L(b)R(X), coefficients N[m][n] of b^m X^n ; divide by (X-b)
    N=[[R[m]*L[n]-L[m]*R[n] for n in range(4)] for m in range(4)]
    # N = (X-b) W, W bidegree (2,2): solve by treating as poly in X with coefficients in b: synthetic division by (X - b)
    # Represent N as poly in X: coeff_n(b)=sum_m N[m][n] b^m
    cX=[[N[m][n] for m in range(4)] for n in range(4)]  # cX[n] is poly in b
    # divide: W_{n-1} = cX[n] + b*W_n  (from top)
    Wn=[None]*3; acc=[Fq(0)]
    for n in range(3,0,-1):
        acc=padd(cX[n], [Fq(0)]+acc)  # cX[n] + b*acc
        Wn[n-1]=acc
    rem=padd(cX[0],[Fq(0)]+acc)
    assert all(r.iszero() for r in rem)
    # Wn[n] poly in b, pad to length 3 (degree<=2)
    out=[]
    for n in range(3):
        w=Wn[n]+[Fq(0)]*5
        assert all(x.iszero() for x in w[3:]), "b-degree >2"
        out.append(w[:3])
    return out  # out[n][m] coefficient of X^n b^m
def resX(W1,W2):
    # homogeneous resultant in X of two binary quadratics (degree 2 in X), coefficients polys in b (deg<=2, homogenized)
    # Sylvester 4x4 with poly entries -> compute by evaluating at 9 points of b and interpolating? use Leibniz expansion.
    a=[W1[n] for n in range(2,-1,-1)]; c=[W2[n] for n in range(2,-1,-1)]
    Z=[Fq(0)]
    S=[[a[0],a[1],a[2],Z],[Z,a[0],a[1],a[2]],[c[0],c[1],c[2],Z],[Z,c[0],c[1],c[2]]]
    tot=[Fq(0)]*9
    for perm in itertools.permutations(range(4)):
        sgn=1
        for i in range(4):
            for j in range(i+1,4):
                if perm[i]>perm[j]: sgn=-sgn
        t=[Fq(1)]
        for i in range(4): t=polmul(t,S[i][perm[i]])
        t=t+[Fq(0)]*(9-len(t))
        tot=[x+(y if sgn>0 else -y) for x,y in zip(tot,t[:9])]
    return tot  # binary octic coefficients low->high in b (degree 8 slot = coefficient of b^8)
def res_bin(f,g):
    # homogeneous resultant of two binary octics given by length-9 coefficient lists
    n=8; F=f[::-1]; G=g[::-1]; Z=Fq(0)
    S=[[Z]*i+F+[Z]*(n-1-i) for i in range(n)]+[[Z]*i+G+[Z]*(n-1-i) for i in range(n)]
    return det(S)
def run(p,ucodes,vcodes):
    Fq.p=p
    u=[dec(c,p) for c in ucodes]; v=[dec(c,p) for c in vcodes]
    ub=[x.conj() for x in u]; vb=[x.conj() for x in v]
    for k in range(1,6):
        su=sum((x**k for x in u),Fq(0)); sv=sum((x**k for x in v),Fq(0))
        assert su.y==sv.y, ("moment",k)
    allr=u+ub+v+vb; assert len({(r.x,r.y) for r in allr})==12
    U=prod(u); V=prod(v); Ub=prod(ub); Vb=prod(vb)
    h=padd(polmul(Ub,V),pneg(polmul(U,Vb)))
    assert all(x.iszero() for x in h[1:]) and not h[0].iszero()
    nz=0; degs=set()
    for pv in itertools.permutations(range(3)):
        for pvb in itertools.permutations(range(3)):
            vv=[v[i] for i in pv]; vvb=[vb[i] for i in pvb]
            Ws=[]
            for j in range(3):
                ks=[k for k in range(3) if k!=j]
                L=prod([vvb[j]]+[u[k] for k in ks]); R=prod([ub[j]]+[vv[k] for k in ks])
                Ws.append(W_coeffs(L,R))
            E12=resX(Ws[0],Ws[1]); E23=resX(Ws[1],Ws[2])
            degs.add((max(i for i,x in enumerate(E12) if not x.iszero()) if any(not x.iszero() for x in E12) else -1,
                      max(i for i,x in enumerate(E23) if not x.iszero()) if any(not x.iszero() for x in E23) else -1))
            r=res_bin(E12,E23)
            if not r.iszero(): nz+=1
    print(f"p={p}: matchings with Res_b(E12,E23)!=0: {nz}/36 ; (deg E12,deg E23) seen {degs}")
run(7,(7,12,26),(8,10,24))
run(11,(11,21,45),(12,13,44))
run(19,(19,233,260),(20,21,106))
