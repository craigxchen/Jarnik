# Part (4): missing rows; cofactor kernel vectors of T=[1;S_M]; exact
from itertools import combinations
from fractions import Fraction
import random, math
from paley1 import chi, S_of
def det(M):
    M=[[Fraction(x) for x in r] for r in M]; n=len(M); d=Fraction(1)
    for c in range(n):
        p=next((r for r in range(c,n) if M[r][c]!=0),None)
        if p is None: return 0
        if p!=c: M[c],M[p]=M[p],M[c]; d=-d
        d*=M[c][c]
        for r in range(c+1,n):
            f=M[r][c]/M[c][c]
            for t in range(c,n): M[r][t]-=f*M[c][t]
    return int(d)
def pivots(T):
    # row-independent set and pivot columns via elimination
    M=[[Fraction(x) for x in r] for r in T]; rows=len(M); cols=len(M[0])
    piv=[]; r=0; order=list(range(rows))
    A=[row[:] for row in M]
    for c in range(cols):
        p=next((i for i in range(r,rows) if A[i][c]!=0),None)
        if p is None: continue
        A[r],A[p]=A[p],A[r]; order[r],order[p]=order[p],order[r]
        for i in range(rows):
            if i!=r and A[i][c]!=0:
                f=A[i][c]/A[r][c]; A[i]=[A[i][t]-f*A[r][t] for t in range(cols)]
        piv.append(c); r+=1
        if r==rows: break
    return order[:r],piv
def kernel_cofactor(T):
    R,piv=pivots(T); r=len(piv); m=len(T[0])
    Tr=[T[i] for i in R]
    M0=[[Tr[a][c] for c in piv] for a in range(r)]
    D=det(M0); out=[]
    for f in range(m):
        if f in piv: continue
        x=[0]*m; x[f]=D
        for j,c in enumerate(piv):
            Mj=[[ (-Tr[a][f] if jj==j else M0[a][jj]) for jj in range(r)] for a in range(r)]
            x[c]=det(Mj)
        out.append(x)
    return r,out
def rank(V):
    if not V: return 0
    _,p=pivots(V); return len(p)
def primes_from(P,cnt):
    out=[];n=P
    while len(out)<cnt:
        if n>1 and all(n%d for d in range(2,int(n**.5)+1)): out.append(n)
        n+=1
    return out
def run(m,kmax,nper,P,seed=1):
    S=S_of(m); random.seed(seed)
    pr=primes_from(P,nper*m); Q=[math.prod(pr[nper*v:nper*v+nper]) for v in range(m)]; q=math.prod(Q)
    # random genuine two-class y with Paley cuts
    def crt(res):
        x=0
        for p,rr in res: 
            Mp=q//p; x+=rr*Mp*pow(Mp,-1,p)
        return x%q
    loc={}
    for v in range(m):
        for p in pr[nper*v:nper*v+nper]:
            r0=random.randrange(p); r1=(r0+random.randrange(1,p))%p; loc[p]=(v,r0,r1)
    y=[crt([(p,(loc[p][1] if S[i][loc[p][0]]==1 else loc[p][2])) for p in pr]) for i in range(m)]
    nch=0; nann=0; bad=0
    for k in range(kmax+1):
        for Mset in combinations(range(m),k):
            nch+=1
            T=[[1]*m]+[S[i] for i in Mset]
            r,X=kernel_cofactor(T)
            assert r==k+1
            imgs=[]
            for x in X:
                nann+=1
                supp=[v for v in range(m) if x[v]!=0]
                assert len(supp)<=r+1 and max(abs(t) for t in x)<=r**(r/2)+1e-9
                a=[sum(S[i][t]*x[t] for t in range(m)) for i in range(m)]
                assert all(a[i]==0 for i in Mset) and sum(a)==0
                for v in range(m):
                    if v not in supp:
                        assert sum(a[i]*S[i][v] for i in range(m))==0
                        assert sum(a[i] for i in range(m) if S[i][v]==1)==0
                    else:
                        assert sum(a[i]*S[i][v] for i in range(m))==(m+1)*x[v]
                div=q//math.prod(Q[v] for v in supp)
                if sum(a[i]*y[i] for i in range(m))%div: bad+=1
                imgs.append(a)
            # images span zero-sum vectors on retained rows: dim m-k-1
            assert rank(imgs)==m-k-1
    return nch,nann,bad
if __name__=="__main__":
    print("m=11,k<=2,1 prime/col",run(11,2,1,101))
    print("m=11,k<=2,2 primes/col",run(11,2,2,101))
    print("m=7,k<=3,2 primes/col",run(7,3,2,53))
    print("m=19,k<=1,1 prime/col",run(19,1,1,211))
