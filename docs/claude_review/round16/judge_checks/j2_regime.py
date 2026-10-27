# Judge check J2: regime arithmetic for sdc_q (n=q+1, d=q, r=2q, q prime = 3 mod 4), 50-digit mpmath.
import mpmath as mp, sympy as sp
mp.mp.dps=50
def par(d,r):
    w=mp.sqrt(8*d*mp.log(4*r)); return w, int(mp.floor(w/2)), int(mp.ceil(4*w))
first=None; first2=None; maxdiff=-10; viol_mlt=0; viol8j4=0; cnt=0
q=3
while q<200000:
    q=sp.nextprime(q)
    if q%4!=3: continue
    d=q; r=2*q
    if d < 316*mp.log(4*r): continue
    cnt+=1
    w,j,m=par(d,r)
    if first is None: first=(q,float(w),j,m)
    if first2 is None and 2*m<q: first2=(q,m)
    if m>=q: viol_mlt+=1
    if m>8*j+4: viol8j4+=1
    maxdiff=max(maxdiff,m-8*j)
print("regime sdc primes q<2e5 (q=3 mod 4):",cnt)
print("first regime sdc q, w, j, m':",first)
print("first q with 2m'<q:",first2)
print("# regime q with m'>=q:",viol_mlt,"; # with m'>8j+4:",viol8j4,"; max m'-8j:",maxdiff)
# Sylvester n=2^k, d=2^k-1, r=d
for k in range(11,16):
    d=2**k-1; r=d; w,j,m=par(d,r)
    print("Syl k=%d d=%d regime=%s w=%.2f j=%d m'=%d"%(k,d,d>=316*mp.log(4*r),float(w),j,m))
