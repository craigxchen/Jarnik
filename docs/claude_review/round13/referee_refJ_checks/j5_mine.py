# independent recomputation of J5 / refJ §1 numbers
import mpmath as mp, math, sys
from fractions import Fraction as Fr
mp.mp.dps=60
def params(k):
    n=2**k; d=n-1
    w=mp.sqrt(8*d*mp.log(4*d))
    j=int(mp.floor(w/2)); mp_=int(mp.ceil(4*w))
    # distances to integers
    dj=min(w/2-mp.floor(w/2), mp.ceil(w/2)-w/2); dm=min(4*w-mp.floor(4*w), mp.ceil(4*w)-4*w)
    return n,d,w,j,mp_,dj,dm
def R1(k,mp_,d,exact=True):
    if exact: return 8*math.comb(mp_+d-1,d) <= 3**d
    return math.log(8)+math.lgamma(mp_+d)-math.lgamma(d+1)-math.lgamma(mp_) <= d*math.log(3)
out=[]
for k in range(12,22):
    n,d,w,j,m,dj,dm=params(k)
    r1 = R1(k,m,d, exact=(k<=20))
    N2=Fr(n*(n-1),2)
    def Q(B):
        s=n//2**B; return N2/(2**B*Fr(s*(s-1),2))
    a = m*mp.log(4)-mp.log(mp.mpf(d)/(m+1))
    o1=m-2*j
    cn=mp.mpf(n)/(n-1)
    def gstar(o,Qv):
        Qv=mp.mpf(Qv.numerator)/Qv.denominator
        if o*Qv-d<=0: return None
        return 2*(a+o*mp.log(cn*Qv))/(o*Qv-d)
    Wmin=3*n*mp.log(n)-14*n; abar_ne=Wmin/(2*d)
    print(f"k={k} n={n} j={j} m'={m} 8j<m':{8*j<m} dist(w/2,Z)={mp.nstr(dj,5)} dist(4w,Z)={mp.nstr(dm,5)} R1={r1} o1={o1} o4={m-8*j}")
    print(f"   Q4={mp.nstr(mp.mpf(Q(4).numerator)/Q(4).denominator,8)} Q5={mp.nstr(mp.mpf(Q(5).numerator)/Q(5).denominator,8)} o1Q4-d={mp.nstr(o1*mp.mpf(Q(4).numerator)/Q(4).denominator-d,8)} o1Q5-d={mp.nstr(o1*mp.mpf(Q(5).numerator)/Q(5).denominator-d,8)} g*(Q4)={gstar(o1,Q(4)) and mp.nstr(gstar(o1,Q(4)),6)} g*(Q5)={gstar(o1,Q(5)) and mp.nstr(gstar(o1,Q(5)),6)} abar_ne={mp.nstr(abar_ne,6)} 2ln15={mp.nstr(2*mp.log(15),6)}")
    # per-subspace with spans (s,c): bits B=2s+c, rho<=(3s+2c)j
    for (s,c) in [(1,0),(2,0),(1,1),(1,2),(0,1),(0,2),(0,3),(0,4),(2,1)]:
        B=2*s+c; o=m-(3*s+2*c)*j
        if o<=0: continue
        Qb=Q(B); Qf=mp.mpf(Qb.numerator)/Qb.denominator
        sc=n//2**B; Pb=2**B*Fr(sc*(sc-1),2)
        if o*Qf-d<=0: continue
        lncoef=mp.log(mp.mpf(n*sc)/ (mp.mpf(Pb.numerator)/Pb.denominator))
        # failure (b): m ln4 + o(lncoef - Qf abar/2) <= ln(d/(m+1)) - d abar/2
        gb=2*(m*mp.log(4)+o*lncoef-mp.log(mp.mpf(d)/(m+1)))/(o*Qf-d)
        delta=mp.mpf(1)/(2**B-1)
        # mixing: (delta-eta)/(1+eta) >= e^{-abar/2}; find smallest abar>=gb satisfying, by bisection on monotone fn
        def ok(ab):
            eta=mp.e**(lncoef-Qf*ab/2)
            if eta>mp.mpf(1)/2: return False
            if not (m*mp.log(4)+o*mp.log(eta) <= mp.log(mp.mpf(d)/(m+1))-d*ab/2): return False
            if (delta-eta)/(1+eta) < mp.e**(-ab/2): return False
            # also P_b <= N2*ab/2 for monotonicity
            if Pb > N2*Fr(int(ab*10**6),10**6)/2: return False
            return True
        lo,hi=mp.mpf(0),mp.mpf(100)
        if not ok(hi): print("   span",(s,c),"no failure <=100"); continue
        for _ in range(200):
            mid=(lo+hi)/2
            if ok(mid): hi=mid
            else: lo=mid
        print(f"   span (s,c)={(s,c)} B={B} o={o} oQ-d={mp.nstr(o*Qf-d,7)} g_b={mp.nstr(gb,6)} exact-threshold~{mp.nstr(hi,7)} W/Wmin={mp.nstr(hi*2*d/Wmin,4)}")
    sys.stdout.flush()
