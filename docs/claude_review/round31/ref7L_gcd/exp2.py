import itertools, math
import numpy as np
box=[w for w in itertools.product(range(-2,3),repeat=3) if any(w)]
pairs=[]
for i,w in enumerate(box):
    for w2 in box[i+1:]:
        c=np.cross(w,w2)
        if c.any(): pairs.append((w,w2))
W1=np.array([p[0] for p in pairs]); W2=np.array([p[1] for p in pairs])
def worst(n):
    a=np.abs(W1@n); b=np.abs(W2@n)
    g=np.gcd(a,b)
    g[(a==0)&(b==0)]=10**12
    return g.max()
print("generic sequence n=(l, round(l*sqrt2), round(l*sqrt3)): max over windows of worst/L")
for lo in [10**2,10**3,10**4,10**5,10**6]:
    best=0
    for l in range(lo,2*lo, max(1,lo//2000)):
        n=np.array([l,round(l*2**.5),round(l*3**.5)]); L=n.max()
        best=max(best,worst(n)/L)
    print(lo,'max ratio in window',round(best,4))
print("adversarial n=(m, m p, m p - m//2 + 1), L=m p, p ~ log2 L : worst/L")
for L in [10**3,10**4,10**5,10**6,10**7,10**8]:
    p=int(math.log2(L)); m=L//p
    n=np.array([m,m*p,m*p-m//2+1])
    print(L,p,round(worst(n)/n.max(),4))
print("hypothesis violated, n=(l,l+1,l+2): F=Z1Z3-Z2^2 vanishes identically -> gcd = deg G(t^n) = L for G=Z1-1")
