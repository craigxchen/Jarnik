import numpy as np, itertools, sys
def chi(a,m):
    a%=m
    if a==0: return 0
    return 1 if pow(a,(m-1)//2,m)==1 else -1
def S(m):
    return np.array([[ (chi(j-i,m) if i!=j else -1) for j in range(m)] for i in range(m)])
for m in [7,11,19,23]:
    s=S(m)
    G=s@s.T; assert (G==(m+1)*np.eye(m,dtype=int)-1).all()
    keep=[i for i in range(m) if s[i,0]==1]; miss=[i for i in range(m) if s[i,0]!=1]
    T=np.vstack([np.ones(m,dtype=int),s[miss]])
    r=np.linalg.matrix_rank(T)
    spark=None
    for t in range(1,r+2):
        combos=list(itertools.combinations(range(m),t))
        found=False
        for c0 in range(0,len(combos),20000):
            cs=np.array(combos[c0:c0+20000])
            sub=T[:,cs].transpose(1,0,2).astype(float)  # (batch, rows, t)
            rk=np.linalg.matrix_rank(sub)
            if (rk<t).any(): found=True; break
        if found: spark=t; break
    print(f"m={m}: retained rows {len(keep)} (density {len(keep)/m:.3f}), T is {T.shape}, rank {r}, spark(T)={spark}, need < m/2={m/2}")
    sys.stdout.flush()
