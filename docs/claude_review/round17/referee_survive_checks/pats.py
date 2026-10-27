# Independent: record sdc_q construction (round5 exact Thm 4.1 / round12 flib conventions) vs survive.md P_a, P'_a; AGL invariance
import sys
def leg(a,q):
    a%=q
    return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
def record_rel_patterns(q):
    # rows: 0 = inf, 1+x = x in F_q ; Paley-I: H[inf,*]=1, H[x,inf]=-1, H[x,y]=[x==y]+chi(y-x)
    idx=['inf']+list(range(q)); M=q+1
    H=[[0]*M for _ in range(M)]
    for a,x in enumerate(idx):
        for b,y in enumerate(idx):
            if x=='inf': H[a][b]=1
            elif y=='inf': H[a][b]=-1
            else: H[a][b]=(1 if x==y else 0)+leg(y-x,q)
    K=H
    if not all(K[i][j]+K[j][i]==(2 if i==j else 0) for i in range(M) for j in range(M)):
        K=[[H[j][i] for j in range(M)] for i in range(M)]
    assert all(K[i][j]+K[j][i]==(2 if i==j else 0) for i in range(M) for j in range(M))
    D=[K[i][0] for i in range(M)]
    K=[[D[i]*K[i][j]*D[j] for j in range(M)] for i in range(M)]
    pats=set()
    for a in range(1,M):
        S=set(x for x in range(M) if K[x][a]==-1)
        Sp=S^{0,a}
        for T in (S,Sp):
            if 0 in T: T=set(range(M))-T
            pats.add(frozenset(t-1 for t in T))  # coordinates x in F_q
    return pats
def survive_patterns(q):
    P=set()
    for a in range(q):
        P.add(frozenset([a]+[x for x in range(q) if leg(a-x,q)==1]))
        P.add(frozenset([a]+[x for x in range(q) if leg(a-x,q)==-1]))
    return P
primes=[p for p in range(7,400) if all(p%d for d in range(2,int(p**.5)+1)) and p%4==3]
for q in primes:
    R=record_rel_patterns(q); S=survive_patterns(q)
    ok=(R==S) and len(R)==2*q
    agl=all(frozenset((s*x+b)%q for x in P) in S for P in S for s in range(1,q) for b in range(q)) if q<=103 else None
    print(q, len(R), 'equal' if ok else 'DIFFER', 'AGL-invariant' if agl else agl, flush=True)
