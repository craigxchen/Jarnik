"""Referee check of Theorem C, Step-1 form (branch independent): for each split prime p, the open arcs
{theta : |sin(theta - arg z)| < (p Norm z)^(-1/2)}, z primitive with Norm z < p, cover R/(pi/2)Z except
exactly the two points arg(pi_p), pi/2 - arg(pi_p). mpmath at 40 digits; touching tolerance 1e-30."""
import sys, json, math
from mpmath import mp, mpf, atan2, asin, sqrt, pi
mp.dps = 40
def split_primes(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [p for p in range(5,n+1) if s[p] and p%4==1]
def two_squares(p):
    for x in range(1,int(p**0.5)+1):
        y2=p-x*x; y=int(round(y2**0.5))
        if y*y==y2: return x,y
TOL = mpf('1e-30'); HP = pi/2
bad = []
for p in split_primes(int(sys.argv[1])):
    arcs = []
    for a in range(1, int(math.isqrt(p))+1):
        for b in range(0, int(math.isqrt(p))+1):
            n = a*a+b*b
            if n >= p or math.gcd(a,b) != 1: continue
            c = atan2(b,a); rho = asin(1/sqrt(mpf(p*n)))
            arcs.append((c-rho, c+rho))
    # unroll circle: shift arcs into [0, HP) start, sweep twice
    ivs = []
    for lo,hi in arcs:
        k = mp.floor(lo/HP); lo -= k*HP; hi -= k*HP
        ivs.append((lo,hi)); ivs.append((lo+HP,hi+HP))
    ivs.sort()
    # sweep over [HP, 2HP)
    cover = ivs[0][1]; touches=[]; gaps=[]
    for lo,hi in ivs[1:]:
        if lo > cover - TOL and cover >= HP and cover < 2*HP:
            if lo - cover > TOL: gaps.append((float(cover-HP), float(lo-HP)))
            else: touches.append(cover - HP)
        if lo > cover + TOL and cover >= HP and cover < 2*HP: pass
        cover = max(cover, hi)
        if lo >= 2*HP: break
    x,y = two_squares(p)
    exp = sorted([atan2(y,x), atan2(x,y)])
    tt = sorted(set([round(float(t),25) for t in touches]))
    ok = (not gaps) and len(tt)==2 and all(abs(mpf(tt[i])-exp[i])<1e-12 for i in range(2))
    if not ok: bad.append(dict(p=p, gaps=gaps[:5], touches=tt, expected=[float(e) for e in exp]))
    print(json.dumps(dict(p=p, narcs=len(arcs), ngaps=len(gaps), touches=tt, ok=ok)), flush=True)
print(json.dumps(dict(summary=True, bad=bad)))
