# Crude random-model count for QS^prod(2, delta, c): P = first s split primes (p = 1 mod 4).
# U(X) = #{P-smooth Gaussian integers (up to units) of norm <= X}, Rankin/saddle-point estimate
#   log U ~ min_sigma [ sigma*log X - 2*sum_p log(1 - p^-sigma) ]   (upper bound, sharp in exponent)
# Expected number of close pairs at norm ~ X: U^2 * X^(-delta)  (|x1-x2| <= |x|^(1-delta) has prob ~ X^(-delta))
import math
def split_primes(s):
    out=[]; n=5
    while len(out)<s:
        if n%4==1 and all(n%q for q in range(3,int(n**0.5)+1,2)): out.append(n)
        n+=4
    return out
def logU(ps, logX):
    best=float('inf')
    for i in range(1,4000):
        sg=i/2000.0
        v=sg*logX-2*sum(math.log(1-p**(-sg)) for p in ps)
        best=min(best,v)
    return best
for s in [50,200,1000,3000]:
    ps=split_primes(s); L=sum(math.log(p) for p in ps)
    for c,d in [(0.25,0.49),(0.49,0.49)]:
        logX=c*L; lu=logU(ps,logX)
        e=2*lu-d*logX
        print(f"s={s:5d} c={c} delta={d}: log prod p={L:9.1f}  log U={lu:8.1f} ({lu/s:.3f} s)  log E[close pairs]={e:9.1f} ({e/s:+.3f} s)")
