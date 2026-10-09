"""Lemma 3.3(3) (sharp.md), large q, level 1: for a random injection y: rows -> Z/(m/2) and parities
sigma, Pr[sum c_x (sigma_x + 2 y_x) = t mod m] <= 4 c_min/m for every target t.  Vectorised MC with
large counts; reports max over trials of (empirical Pr)/(4 c_min/m) and a 4-sigma upper envelope."""
import math, random
import numpy as np
rng = random.Random(9); g = np.random.default_rng(9)
def primes(n):
    s = bytearray([1])*(n+1); s[0:2]=b"\x00\x00"
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
P=[p for p in primes(400) if p>60]
worst=0; worst_z=-9
for trial in range(120):
    q=rng.choice(P); m=q-1 if q%4==1 else q+1; h=m//2
    R=rng.randint(4, h//2)          # rows, so that m/2 >= 2R
    sigma=np.array([rng.randint(0,1) for _ in range(R)])
    cmax=rng.choice([1,1,2,3,4,6])
    c=np.array([rng.randint(-cmax,cmax) for _ in range(R)]); c[0]-=c.sum()
    if not c.any(): continue
    cmin=int(np.min(np.abs(c[c!=0])))
    N=200000
    # random injections: argsort of random keys, take first R
    keys=g.random((N,h))
    y=np.argpartition(keys, R, axis=1)[:,:R]
    kap=((sigma+2*y)*c).sum(axis=1)%m
    cnt=np.bincount(kap, minlength=m)
    p_emp=cnt.max()/N; bound=4*cmin/m
    z=(cnt.max()-N*bound)/math.sqrt(N*bound)
    worst=max(worst,p_emp/bound); worst_z=max(worst_z,z)
print(f"max empirical/bound over 120 trials (N=2e5 each, all targets t): {worst:.3f};  max z-score vs bound: {worst_z:.2f}")
