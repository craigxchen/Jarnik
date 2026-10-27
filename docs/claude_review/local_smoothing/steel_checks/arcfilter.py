# S1 check: the paper's imaginary light-ray translation, restricted to one circle |xi|=R,
# multiplies mode xi by exp(-l*(|xi| - xi.omega)) = exp(-l R (1-cos(theta-omega))).
# Compare the filtered sum S_omega (l=1) with exact arc counts (arc length C sqrt R), and
# check the Bessel/Hecke expansion S_omega = sum_k I_k(l R) e^{-l R} e^{-ik omega} sum_xi e^{ik theta}.
import math, cmath, itertools
from mpmath import mp, besseli, mpf
mp.dps = 30
def pts(N):
    out=[]; a=0
    while a*a<=N:
        b2=N-a*a; b=int(math.isqrt(b2))
        if b*b==b2:
            for s1 in (1,-1):
                for s2 in (1,-1):
                    out.append((s1*a,s2*b))
        a+=1
    return sorted(set(out))
def run(N,C=1.0,l=1.0):
    P=pts(N); R=math.sqrt(N); th=sorted(math.atan2(y,x)%(2*math.pi) for x,y in P)
    w=C/math.sqrt(R)  # angular length of arc of length C sqrt R
    # exact max arc count
    ext=th+[t+2*math.pi for t in th]; M=0; j=0
    for i in range(len(th)):
        while ext[j]<=ext[i]+w: j+=1
        M=max(M,j-i)
    # filtered max over omega at cluster points (enough: sample omega at points and midpoints)
    best=0; bestw=None
    cands=th+[(th[i]+th[(i+1)%len(th)])/2 for i in range(len(th))]
    for om in cands:
        S=sum(math.exp(-l*R*(1-math.cos(t-om))) for t in th)
        if S>best: best,bestw=S,om
    # Bessel/Hecke expansion at bestw
    lR=mpf(l*R); K=int(6*math.sqrt(l*R))+8
    ser=mpf(0)
    for k in range(-K,K+1):
        if k%4: continue
        h=sum(cmath.exp(1j*k*t) for t in th)
        ser+= besseli(abs(k),lR)*mp.e**(-lR)*mpf((h*cmath.exp(-1j*k*bestw)).real)
    return len(P),M,best,float(ser)
for N in [5**2*13*17*29*37*41, 5*13*17*29*37*41*53*61, 65**4*17, 1105**3]:
    print(N, run(N))
