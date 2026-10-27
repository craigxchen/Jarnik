# Independent construction of sdc_q patterns; symmetry and slack checks. Exact integer arithmetic.
import sys, itertools, numpy as np
def chi(a,q):
    a%=q
    return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
def sdc(q):
    # rows: 0 = inf (base), 1+x for x in F_q. relative patterns (omit base)
    P=[]; Pp=[]
    for a in range(q):
        P.append(frozenset([a]+[x for x in range(q) if chi(a-x,q)==1]))
        Pp.append(frozenset([a]+[x for x in range(q) if chi(a-x,q)==-1]))
    return P,Pp
def check(q):
    P,Pp=sdc(q); pats=set(P)|set(Pp)
    assert len(pats)==2*q
    # translations and multiplications
    trans=all(frozenset((x+1)%q for x in S) in pats for S in pats)
    mult={t:all(frozenset((t*x)%q for x in S) in pats for S in pats) for t in range(1,q)}
    allmult=all(mult.values())
    # pattern incidence with base column: full rows n=q+1, base index q
    n=q+1
    M=np.zeros((2*q,n),dtype=np.int64)
    for i,S in enumerate(sorted(pats,key=sorted)):
        for x in S: M[i,x]=1
    # sigma for pairs incl base: lambda = e_x - e_y
    sig=set()
    for x in range(n):
        for y in range(x+1,n):
            v=M[:,x]-M[:,y]
            sig.add(int(np.abs(v).sum())*2-2*q)  # 2*sigma = sum(2|l(S)|-1)
    return trans,allmult,sorted(s/2 for s in sig),M
if __name__=="__main__":
    for q in [7,11,19,23,31,43,47,59,67,71,79,83]:
        tr,mu,sg,_=check(q)
        print(f"q={q} translation-invariant={tr} all F_q^* multiplications preserve={mu} pair sigmas={sg}",flush=True)
