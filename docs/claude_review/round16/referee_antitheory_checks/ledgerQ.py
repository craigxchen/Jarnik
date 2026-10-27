# crossover n (r=n, C=1, W=Wmin) where L0 - Lc - m(Q abar/2 + ln2) > 0, exponent m, for Q-class strata
import mpmath as mp
mp.mp.dps=30
def margin(n,Q):
    r=n; d=n-1; w=mp.sqrt(8*d*mp.log(4*r)); m=int(mp.ceil(4*w))
    W=3*n*mp.log(n)-14*n; abar=W/(2*d)
    L0=W/4-mp.log(mp.mpf(d)/(m+1)); Lc=m*mp.log(2*d*mp.e/m)
    return L0-Lc-m*(Q*abar/2+mp.log(2))
for Q in [16,32,64,128,256]:
    lo,hi=2900,10**15
    assert margin(hi,Q)>0
    if margin(lo,Q)>0: print(Q,'<=2900'); continue
    while hi-lo>max(1,lo//1000):
        mid=(lo+hi)//2
        if margin(mid,Q)>0: hi=mid
        else: lo=mid
    print('Q=',Q,'crossover n ~',hi, ' margin at 1e6:', mp.nstr(margin(10**6,Q),5))
