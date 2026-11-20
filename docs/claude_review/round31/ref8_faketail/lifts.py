# Own check of the exact integer moment lifts d^(s) (Part VIII.3), exact integers.
import random, math, sys
def paley(q):
    sq={x*x%q for x in range(1,q)}
    chi=[0]+[1 if x in sq else -1 for x in range(1,q)]
    return [[chi[(j-i)%q] for j in range(q)] for i in range(q)]
def solve_mod(M,b,p):
    # solve M x = b mod p (M list of rows), return one solution or None
    n=len(M); m=len(M[0])
    A=[row[:]+[bb] for row,bb in zip(M,b)]
    piv=[];r=0
    for c in range(m):
        pr=next((k for k in range(r,n) if A[k][c]%p),None)
        if pr is None: continue
        A[r],A[pr]=A[pr],A[r]
        inv=pow(A[r][c],-1,p); A[r]=[x*inv%p for x in A[r]]
        for k in range(n):
            if k!=r and A[k][c]%p:
                f=A[k][c]; A[k]=[(x-f*y)%p for x,y in zip(A[k],A[r])]
        piv.append(c); r+=1
        if r==n: break
    for k in range(r,n):
        if A[k][-1]%p: return None
    x=[0]*m
    for k,c in enumerate(piv): x[c]=A[k][-1]
    return x
def rnd(r,q):  # mean-zero rounding of residue r in [0,q)
    return r if random.random()<(q-r)/q else r-q
def run(q,rho,seed):
    random.seed(seed)
    K=paley(q); Acls=[i for i in range(q) if K[i][0]==1]
    assert len(Acls)==(q-1)//2
    S=[i for i in Acls if random.random()<rho]; Sc=[i for i in range(q) if i not in S]
    # identities
    for i in range(q):
        for k in range(q):
            assert sum(K[i][j]*K[k][j] for j in range(q))==(q-1 if i==k else -1)
    J=q//4; worst=0; worstt=0
    for s in range(1,J+1):
        best=None
        for trial in range(20):
            d=[rnd(pow(j,s,q),q) for j in range(q)]
            Kd=[sum(K[i][j]*d[j] for j in range(q)) for i in range(q)]
            assert all(x%q==0 for x in Kd), "tau not integral"
            tau=[x//q for x in Kd]
            mt=max(abs(tau[i]) for i in range(q))
            if best is None or mt<best[0]: best=(mt,d,tau)
        mt,d,tau=best
        rhs=[(-sum(K[i][j]*tau[i] for i in S))%q for j in range(q)]
        M=[[K[i][j]%q for i in Sc] for j in range(q)]
        x=solve_mod(M,rhs,q); assert x is not None, "completion impossible"
        bestd=None
        for trial in range(20):
            tp=[0]*q
            for i in S: tp[i]=tau[i]
            for c,i in enumerate(Sc): tp[i]=rnd(x[c],q)
            KTt=[sum(K[i][j]*tp[i] for i in range(q)) for j in range(q)]
            assert all(v%q==0 for v in KTt)
            ds=[d[j]-KTt[j] for j in range(q)]
            n=max(map(abs,ds))
            if bestd is None or n<bestd[0]: bestd=(n,ds)
        n,ds=bestd
        assert all((ds[j]-pow(j,s,q))%q==0 for j in range(q))
        img=[sum(K[i][j]*ds[j] for j in range(q)) for i in S]
        assert len(set(img))==1
        # also constant after deleting column 0
        img2=[sum(K[i][j]*ds[j] for j in range(1,q)) for i in S]
        assert len(set(img2))==1
        worst=max(worst,n); worstt=max(worstt,mt)
    # s = (q-1)/2 must fail integrality (degree q-1): negative control
    s=(q-1)//2; d=[rnd(pow(j,s,q),q) for j in range(q)]
    Kd=[sum(K[i][j]*d[j] for j in range(q)) for i in range(q)]
    neg=any(x%q for x in Kd)
    print(f"q={q} |S|={len(S)} J={J}: all lifts OK; max|tau|={worstt} "
          f"(sqrt(q log q)={math.sqrt(q*math.log(q)):.1f}); max|d^(s)|={worst} "
          f"(q^1.5 sqrt(log q)={q**1.5*math.sqrt(math.log(q)):.0f}); ratio={worst/(q**1.5*math.sqrt(math.log(q))):.3f}; "
          f"neg control s=(q-1)/2 nonintegral: {neg}")
for q in (103,127):
    run(q,0.99,1)
