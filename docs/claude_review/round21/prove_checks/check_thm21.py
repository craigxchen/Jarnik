"""Exact check of Theorem 2.1 (items 1,2,4) and window statistics on record clusters.
For each n and each window of k>=5 consecutive points with chord constant C<=CMAX: normalise by the Gaussian gcd,
compute N', s=omega(N'), Y=P_+(N'), pi(Y;4,1), W', A=W'/(k log k), and check
  (1) W' < s*(4W'/k + 8 log^+ C)        (2) k < 4 pi(Y;4,1)(1 + 2k log^+C/W')
Also prints sum_j p_j^(-1/2) over primes of N' and whether Y <= k^2 (window at delta=0)."""
import sys, math
from sympy import factorint, primerange
from glib import windows, ggcd, gnorm
KNOWN=[235370525,537606725,1176852625,2163838625,10414625,412333805,1762645]
def pi41(y): return sum(1 for p in primerange(2,y+1) if p%4==1)
CMAX=float(sys.argv[1]) if len(sys.argv)>1 else 9.0
viol=0; tested=0; best={}
for n in KNOWN:
    res,f=windows(n,5,CMAX)
    for k,C,W in res:
        g=W[0]
        for z in W[1:]: g=ggcd(g,z)
        Nn=n//gnorm(g); fn=factorint(Nn)
        Wp=math.log(Nn); s=len(fn); Y=max(fn); C=float(C)
        lc=max(0.0,math.log(C))
        ok1 = Wp < s*(4*Wp/k+8*lc)
        ok2 = k < 4*pi41(Y)*(1+2*k*lc/Wp)
        tested+=1; viol+= (not ok1) + (not ok2)
        key=(n,k)
        if key not in best or C<best[key][0]:
            best[key]=(C,Nn,s,Y,pi41(Y),Wp/(k*math.log(k)),sum(p**-0.5 for p in fn), Y<=k*k)
print("windows tested:",tested,"violations:",viol)
print("n k minC N' omega P+ pi(P+;4,1) A=W'/(k log k) sum p^-1/2 window(P+<=k^2)")
for (n,k),v in sorted(best.items()):
    print(n,k,"%.6f"%v[0],v[1],v[2],v[3],v[4],"%.3f"%v[5],"%.3f"%v[6],v[7])
