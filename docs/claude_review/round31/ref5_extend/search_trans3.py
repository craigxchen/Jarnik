# Char 3: points z_S = prod_{w in S} (Y + w), Y = i*T (conj Y = -Y), w in F_3[t]; S a transversal of a union W of pairs {w,-w}
# (w != 0), up to `maxdrop` pairs of V\{0} dropped.  conj z_S = +-z_{-S}, so all transversals of W have the same norm.
# Group by (W, top coefficients above floor(D/2)) after making monic.  Report least D for each cluster size m>=3.
import sys, itertools, random
# F_9 = F_3[i], element (a,b) -> a+b i, encoded as a+3b (0..8)
def f9mul(x,y):
    a,b=x%3,x//3; c,d=y%3,y//3
    return ((a*c-b*d)%3)+3*((a*d+b*c)%3)
def f9add(x,y): return (x%3+y%3)%3+3*((x//3+y//3)%3)
MUL=[[f9mul(x,y) for y in range(9)] for x in range(9)]
ADD=[[f9add(x,y) for y in range(9)] for x in range(9)]
INV=[None]+[next(y for y in range(1,9) if MUL[x][y]==1) for x in range(1,9)]
def pmul(p,q):
    r=[0]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        if x:
            row=MUL[x]
            for j,y in enumerate(q):
                if y: r[i+j]=ADD[r[i+j]][row[y]]
    return r
def span3(basis):  # basis: list of F_3[t] coefficient lists
    L=max(len(b) for b in basis); S=[tuple([0]*L)]
    for b in basis:
        bb=list(b)+[0]*(L-len(b))
        S=[tuple((s[k]+c*bb[k])%3 for k in range(L)) for s in S for c in range(3)]
    return S
def search(V,T,best,tag,maxdrop):
    L=max(len(T),max(len(v) for v in V))
    nz=[v for v in V if any(v)]
    reps=[];seen=set()
    for v in nz:
        if v in seen: continue
        neg=tuple((-c)%3 for c in v); seen|={v,neg}; reps.append((v,neg))
    def factor(v):  # Y + v = i*T + v, as F_9 coefficient list
        f=[0]*L
        for k,c in enumerate(T): f[k]=ADD[f[k]][MUL[3][c]]  # 3 encodes i
        for k,c in enumerate(v): f[k]=ADD[f[k]][c]
        while len(f)>1 and f[-1]==0: f.pop()
        return f
    fac=[(factor(a),factor(b)) for a,b in reps]
    n=len(reps); groups={}
    def dfs(k,mask,drops,prod):
        if k==n:
            if mask:
                lc=INV[prod[-1]]; p=[MUL[lc][c] for c in prod]; D=len(p)-1; K=D//2
                key=(mask,tuple(p[K+1:]))
                groups.setdefault(key,set()).add(tuple(p))
            return
        if drops<maxdrop: dfs(k+1,mask,drops+1,prod)
        dfs(k+1,mask|1<<k,drops,pmul(prod,fac[k][0]))
        dfs(k+1,mask|1<<k,drops,pmul(prod,fac[k][1]))
    dfs(0,0,0,[1])
    for key,ps in groups.items():
        m=len(ps)
        if m<3: continue
        # recompute D from mask: sum of factor degrees
        D=sum(len(fac[k][0])-1 for k in range(n) if key[0]>>k&1)
        if m not in best or D<best[m][0]: best[m]=(D,tag,bin(key[0]))
    return best
if __name__=="__main__":
    mode=sys.argv[1]; best={}
    if mode=="full2":   # V=F_3[t]_{<2}, all pairs droppable
        V=span3([[1],[0,1]])
        Ts=[list(x) for d in range(0,5) for x in itertools.product(range(3),repeat=d)]+[[0]*L+[1] for L in range(2,13)]
        for T in Ts:
            if not any(T): continue
            search(V,T,best,("F3[t]<2","T="+str(T)),4)
    elif mode=="full3":
        V=span3([[1],[0,1],[0,0,1]])
        Ts=[[0]*L+[1] for L in range(3,int(sys.argv[3]))]+[list(x) for x in itertools.product(range(3),repeat=3)]
        for T in Ts:
            if not any(T): continue
            search(V,T,best,("F3[t]<3","T="+str(T)),int(sys.argv[2]))
    elif mode=="rand":
        random.seed(int(sys.argv[2]))
        for it in range(int(sys.argv[3])):
            n=random.randint(2,3); e=random.randint(n,5)
            basis=[[1]]+[[random.randrange(3) for _ in range(e)] for _ in range(n-1)]
            V=span3(basis)
            if len(set(V))<3**n: continue
            T=random.choice([[random.randrange(3) for _ in range(random.randint(1,6))],[0]*random.randint(e,10)+[1]])
            search(V,T,best,("basis="+str(basis),"T="+str(T)),2 if n==3 else 4)
    for m in sorted(best): print(m,best[m])
