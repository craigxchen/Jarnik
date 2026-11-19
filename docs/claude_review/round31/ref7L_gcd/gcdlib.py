import numpy as np, random, math
P = 2147483629  # prime < 2^31
def trim(a):
    i=len(a)
    while i>0 and a[i-1]==0: i-=1
    return a[:i]
def pgcd(a,b):
    """gcd degree of polys (numpy int64 arrays, low->high coeffs) mod P"""
    a=trim(np.array(a,dtype=np.int64)%P); b=trim(np.array(b,dtype=np.int64)%P)
    if len(a)<len(b): a,b=b,a
    while len(b)>0:
        # a mod b
        inv=pow(int(b[-1]),P-2,P)
        db=len(b)-1
        a=a.copy()
        for i in range(len(a)-1,db-1,-1):
            c=int(a[i])
            if c:
                f=(c*inv)%P
                a[i-db:i+1]=(a[i-db:i+1]-f*b)%P
        a=trim(a[:db])
        a,b=b,a
    return len(a)-1
def spec(F,n):
    """F: dict exponent-tuple->coef. returns dense poly in t after clearing min power"""
    e={}
    for u,c in F.items():
        k=sum(x*y for x,y in zip(u,n)); e[k]=(e.get(k,0)+c)%P
    e={k:v for k,v in e.items() if v}
    if not e: return None
    m=min(e); M=max(e); a=np.zeros(M-m+1,dtype=np.int64)
    for k,v in e.items(): a[k-m]=v
    return a
def common(Fs,n):
    polys=[spec(F,n) for F in Fs]
    polys=[p for p in polys if p is not None]
    g=None
    for p in polys:
        if g is None: g=p
        else:
            d=pgcd(g,p); 
            # recompute actual gcd poly needed for k>2; we only use k=2
            return d
    return len(g)-1
