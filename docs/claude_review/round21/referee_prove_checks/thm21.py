# Referee check of Thm 2.1(1)-(2) consequences on all windows (M>=5 consecutive points, chord constant <= 9)
import math, sys
sys.argv=['x']
exec(open(__file__.replace('thm21.py','records.py')).read().split('for (name,n,k) in')[0])
from sympy import factorint, primerange
def pi41(Y): return sum(1 for p in primerange(2,Y+1) if p%4==1)
tot=0; viol=0
for n in [1176852625,235370525,537606725,2163838625]:
    P=points(n); m=len(P); R=math.sqrt(n)
    for M in range(5,m):
        for i in range(m):
            win=[P[(i+t)%m] for t in range(M)]
            ch=math.sqrt(gnorm((win[-1][0]-win[0][0],win[-1][1]-win[0][1])))
            # chord constant: max pairwise chord / sqrt R  (arc endpoints give the max for an arc < pi)
            a0=math.atan2(win[0][1],win[0][0]); a1=math.atan2(win[-1][1],win[-1][0]); span=(a1-a0)%(2*math.pi)
            if span>=math.pi: continue
            C=ch/math.sqrt(R)
            if C>9: continue
            G=win[0]
            for z in win[1:]: G=ggcd(G,z)
            Np=n//gnorm(G); f=factorint(Np); s=len(f); Y=max(f); Wp=math.log(Np)
            lc=max(0.0,math.log(C))
            ok1 = Wp < s*(4*Wp/M+8*lc)
            ok2 = M < 4*pi41(Y)*(1+2*M*lc/Wp)
            tot+=1; viol+= (not ok1) or (not ok2)
            if (not ok1) or (not ok2): print('VIOL',n,M,i,C,f)
print("windows checked:",tot,"violations:",viol)
