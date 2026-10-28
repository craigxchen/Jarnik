# Referee's own modular-rank reg; prime p=2147483647 (different from width.md's 2^61-1).
import numpy as np, itertools, sys
P=2147483647
def T(n,a,b):
    out=[]
    def rec(pre,k,s,pos):
        if k==n-1:
            y=a-b-s; pos2=pos+max(y,0)
            if pos2<=a: out.append(tuple(pre+[y]))
            return
        for c in range(-b,a+1):
            if pos+max(c,0)<=a: rec(pre+[c],k+1,s+c,pos+max(c,0))
    rec([],0,0,0); return out
def rank_mod(A):
    A=A.copy()%P; r=0; rows,cols=A.shape
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if A[i,c]: piv=i;break
        if piv is None: continue
        A[[r,piv]]=A[[piv,r]]
        inv=pow(int(A[r,c]),P-2,P)
        A[r]=(A[r]*inv)%P
        nz=np.nonzero(A[:,c])[0]
        for i in nz:
            if i!=r: A[i]=(A[i]-A[i,c]*A[r])%P
        r+=1
        if r==rows: break
    return r
def exps(k,d):
    for tot in range(d+1):
        for c in itertools.combinations_with_replacement(range(k),tot):
            e=[0]*k
            for j in c:e[j]+=1
            yield e
def full_at(S,d):
    pts=np.array([s[:-1] for s in S],dtype=object)  # drop last coord (affine iso on hyperplane)
    k=pts.shape[1]
    if k==0: return len(S)<=1
    cols=[]
    for e in exps(k,d):
        col=[1]*len(S)
        for i,s in enumerate(S):
            v=1
            for j in range(k): v=v*pow(s[j]%P,e[j],P)%P
            col[i]=v
        cols.append(col)
    A=np.array(cols,dtype=np.int64).T
    return rank_mod(A)==len(S)
def regp(S,start=0):
    d=start
    while not full_at(S,d): d+=1
    return d
if __name__=="__main__":
    pass
