# Max over heights X >= (prod p)^c of the random-model count log E(X) = 2 log U(X) - delta log X,
# P = first s split primes; Rankin/saddle estimate for log U. Also report the saddle sigma*.
import math, sys
def split_primes(s):
    out=[]; n=5
    while len(out)<s:
        if n%4==1 and all(n%q for q in range(3,int(n**0.5)+1,2)): out.append(n)
        n+=4
    return out
def logU(lps, logX):
    best=(float('inf'),None)
    lo,hi=1e-4,3.0
    f=lambda sg: sg*logX-2*sum(math.log(1-math.exp(-sg*l)) for l in lps)
    for _ in range(100):  # ternary search, f convex in sigma
        m1=lo+(hi-lo)/3; m2=hi-(hi-lo)/3
        if f(m1)<f(m2): hi=m2
        else: lo=m1
    sg=(lo+hi)/2; return f(sg), sg
for s in [50,200,1000,3000,10000]:
    lps=[math.log(p) for p in split_primes(s)]; L=sum(lps)
    for c,d in [(0.25,0.49),(0.49,0.49)]:
        best=None
        for k in range(0,200):
            t=c*(1.05**k)            # log X = t*L
            lu,sg=logU(lps,t*L); e=2*lu-d*t*L
            if best is None or e>best[0]: best=(e,t,sg,lu)
        e,t,sg,lu=best
        print(f"s={s:6d} c={c} delta={d}: max log E = {e:10.1f} ({e/s:+.3f} s) at log X = {t:.2f} log prod p, sigma*={sg:.3f}")
    sys.stdout.flush()
