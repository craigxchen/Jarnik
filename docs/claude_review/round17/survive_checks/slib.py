# survive_checks library: sdc_q box, symmetry, monomial bases. Independent code (round 17, survival attempt).
import itertools, numpy as np
def chi(a,q):
    a%=q
    return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
def patterns(q):
    """relative sdc_q patterns (base row omitted): P_a={a}u{x:chi(a-x)=1}, P'_a={a}u{x:chi(a-x)=-1}"""
    P=[frozenset([a]+[x for x in range(q) if chi(a-x,q)==1]) for a in range(q)]
    Pp=[frozenset([a]+[x for x in range(q) if chi(a-x,q)==-1]) for a in range(q)]
    return P+Pp
def box(q,L):
    """all k in {-1,0,1}^q with |k(P)|<=L for all patterns; enumerated by DFS with pruning"""
    pats=patterns(q)
    inc=np.zeros((q,2*q),dtype=np.int64)
    for i,S in enumerate(pats):
        for x in S: inc[x,i]=1
    # remaining-count per pattern after position x
    rem=np.zeros((q+1,2*q),dtype=np.int64)
    for x in range(q-1,-1,-1): rem[x]=rem[x+1]+inc[x]
    out=[]; k=[0]*q
    def dfs(x,s):
        if x==q:
            out.append(tuple(k)); return
        for v in (0,1,-1):
            s2=s+v*inc[x]
            # feasibility: need |s2|-rem[x+1] <= L  (can still be corrected)
            if np.all(np.abs(s2)-rem[x+1]<=L) :
                k[x]=v; dfs(x+1,s2)
        k[x]=0
    dfs(0,np.zeros(2*q,dtype=np.int64))
    A=[a for a in out if all(abs(sum(a[x] for x in S))<=L for S in pats)]
    return A
if __name__=="__main__":
    import time
    for q in [7,11,19,23]:
        for L in [1,2,3]:
            t=time.time()
            if q>=19 and L>=3: continue
            A=box(q,L); print(q,L,len(A),round(time.time()-t,1),flush=True)
