"""Appendix A (sharp.md) constants.
(1) Lemma A.2: ||1_(f=v)||_A <= prod_(u != v, |u|<=3) (3+|u|)/|v-u| for v in {-3..3}; max over v != 0.
(2) Verify the Wiener-norm level-set bound on random integer f: F_2^m -> Z with ||f||_A <= 3
    (built as signed sums of few subgroup/coset indicators), m = 4..6: ||1_(f=v)||_A <= bound(v).
(3) Gaussian-binomial constant prod_(i>=1) (1 - 2^-i)^-1 used in '#affine subspaces <= 3.47 (m+1) 2^(m^2/4+m)'.
(4) Lemma A.1 'support 4 aligned' probability: Pr[pi(x4) = pi(x1)+pi(x2)+pi(x3)] = 1/(M-3), by MC.
"""
import itertools, random, math
import numpy as np
rng = random.Random(4)
out=[]
bounds={}
for v in range(-3,4):
    b=1.0
    for u in range(-3,4):
        if u!=v: b*= (3+abs(u))/abs(v-u)
    bounds[v]=b
out.append("(1) level-set bounds: "+", ".join(f"v={v}: {bounds[v]:.0f}" for v in bounds)+f";  max over v != 0 = {max(bounds[v] for v in bounds if v):.0f}")
def walsh(m):
    n=1<<m
    H=np.array([[(-1)**bin(a&x).count('1') for x in range(n)] for a in range(n)])
    return H
worst=0; cnt=0
for m in (4,5,6):
    n=1<<m; H=walsh(m)
    # random subspaces/cosets
    def rand_coset():
        d=rng.randint(0,m)
        basis=[rng.randrange(1,n) for _ in range(d)]
        span={0}
        for bvec in basis:
            span|={s^bvec for s in span}
        x0=rng.randrange(n)
        ind=np.zeros(n,dtype=np.int64)
        for s in span: ind[s^x0]=1
        return ind
    for it in range(4000):
        f=np.zeros(n,dtype=np.int64)
        for _ in range(rng.randint(1,4)):
            f+= rng.choice([1,-1])*rand_coset()
        fh=H@f/n
        A=np.abs(fh).sum()
        if A>3+1e-9: continue
        cnt+=1
        for v in set(f.tolist()):
            if v==0: continue
            ind=(f==v).astype(np.int64)
            Av=np.abs(H@ind/n).sum()
            worst=max(worst, Av/bounds[v])
out.append(f"(2) {cnt} random integer f with ||f||_A <= 3 (m=4..6): max ||1_(f=v)||_A / bound(v) = {worst:.4f} (<= 1 required); values lie in [-3,3]")
c=1.0
for i in range(1,200): c/= (1-2.0**-i)
out.append(f"(3) prod (1-2^-i)^-1 = {c:.5f}")
for m in (5,7):
    M=1<<m; N=200000; hits=0
    for _ in range(N):
        a,b2,c3,d=rng.sample(range(M),4)
        hits += (a^b2^c3)==d
    out.append(f"(4) M={M}: MC Pr[4 random distinct points coplanar] = {hits/N:.5f} vs 1/(M-3) = {1/(M-3):.5f}")
print("\n".join(out))
