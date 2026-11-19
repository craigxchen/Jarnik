"""Own F_9[t] gcd (F_9 = F_3[i], code a+3b <-> a+bi), table-driven numpy Euclid.
Used on the ternary referee's saved four-mark reduced points (data only)."""
import numpy as np, itertools, sys
ADD=np.zeros((9,9),np.int64); MUL=np.zeros((9,9),np.int64)
for x in range(9):
    for y in range(9):
        a,b=x%3,x//3; c,d=y%3,y//3
        ADD[x,y]=((a+c)%3)+3*((b+d)%3)
        MUL[x,y]=((a*c-b*d)%3)+3*((a*d+b*c)%3)
NEG=np.array([((-(x%3))%3)+3*((-(x//3))%3) for x in range(9)])
INV=np.zeros(9,np.int64)
for x in range(1,9):
    for y in range(1,9):
        if MUL[x,y]==1: INV[x]=y
def trim(a):
    nz=np.nonzero(a)[0]; return a[:nz[-1]+1] if len(nz) else a[:0]
def polymod(a,b):
    a=a.copy(); db=len(b)-1; ib=INV[b[-1]]
    nb=NEG[b]
    for k in range(len(a)-1,db-1,-1):
        c=a[k]
        if c:
            f=MUL[c,ib]
            a[k-db:k+1]=ADD[a[k-db:k+1],MUL[f,nb]]
    return trim(a[:db])
def gcd(a,b):
    a=trim(a); b=trim(b)
    while len(b):
        a,b=b,polymod(a,b)
    return a
def load(path):
    arr=np.load(path)
    pts=[]
    for row in arr:
        h=len(row)//2; re,im=row[:h]%3,row[h:]%3
        pts.append(trim(re+3*im))
    return pts
if __name__=='__main__':
    pts=load(sys.argv[1]); m=len(pts); D=len(pts[0])-1
    print("m",m,"D",D,[len(p)-1 for p in pts])
    g={}
    for r in range(2,m+1):
        for T in itertools.combinations(range(m),r):
            x=pts[T[0]]
            for j in T[1:]: x=gcd(x,pts[j])
            g[T]=len(x)-1
            print("gcd",T,g[T],flush=True)
    for a,b in itertools.combinations(range(m),2):
        d=trim(np.array([ADD[x,NEG[y]] for x,y in zip(pts[a],pts[b])]))
        print("diff",(a,b),len(d)-1,"chord deg c_ab =",len(d)-1-g[(a,b)])
