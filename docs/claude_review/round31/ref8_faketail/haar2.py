# Fast exact version of the tangent-space check: containment checked over Z,
# ranks mod a large prime p give rank_Q(V) >= r_p(V) and dimT <= |I|-r_p(Dm); equality closes it.
import random
p=(1<<61)-1
def rank_mod(M):
    M=[r[:] for r in M]; r=0; n=len(M); m=len(M[0]) if M else 0
    for c in range(m):
        pr=next((k for k in range(r,n) if M[k][c]%p),None)
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]; inv=pow(M[r][c],-1,p)
        M[r]=[x*inv%p for x in M[r]]
        for k in range(r+1,n):
            if M[k][c]%p:
                f=M[k][c]; M[k]=[(x-f*y)%p for x,y in zip(M[k],M[r])]
        r+=1
    return r
def paley(q):
    sq={x*x%q for x in range(1,q)}
    chi=[0]+[1 if x in sq else -1 for x in range(1,q)]
    return [[chi[(j-i)%q] for j in range(q)] for i in range(q)]
for q,rho,seed in ((43,0.9,1),(59,0.95,2),(103,0.99,3),(127,0.99,4)):
    random.seed(seed)
    K=paley(q); A=[i for i in range(q) if K[i][0]==1]
    S=[i for i in A if random.random()<rho]; Sc=[i for i in range(q) if i not in S]
    I=list(range(1,q))
    Dm=[[K[i][j]-K[S[0]][j] for j in I] for i in S[1:]]
    Vcols=[[1]*len(I)]+[[K[i][j] for j in I] for i in Sc]
    contained=all(sum(r[c]*v[c] for c in range(len(I)))==0 for r in Dm for v in Vcols)
    rD=rank_mod(Dm); rV=rank_mod(Vcols)
    # normalisation: delete columns constant on S, project
    keep=[c for c,j in enumerate(I) if len({K[i][j] for i in S})>1]
    Dm2=[[r[c] for c in keep] for r in Dm]; V2=[[v[c] for c in keep] for v in Vcols]
    rD2=rank_mod(Dm2); rV2=rank_mod(V2)
    print(f"q={q} |S|={len(S)} contained={contained} rank(Dm)={rD} dimT<= {len(I)-rD} rank V={rV} "
          f"EQUAL={contained and rV==len(I)-rD}; after deleting {len(I)-len(keep)} constant cols: "
          f"projection onto T' surjective={rV2==len(keep)-rD2}")
