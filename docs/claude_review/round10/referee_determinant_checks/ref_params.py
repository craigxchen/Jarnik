# Referee: independent exact check of Theorem 2.1 parameters and Theorem 4.1 constants.
# m = largest integer with 8*binom(m+n-2, n-1) <= (2T+1)^(n-1), T = n  (exact big-integer arithmetic).
import math, sys
from math import comb, log, sqrt

def exact_m(n, T):
    target = (2*T+1)**(n-1)
    ok = lambda m: 8*comb(m+n-2, n-1) <= target
    # start from lgamma estimate
    def lb(a,b): return math.lgamma(a+1)-math.lgamma(b+1)-math.lgamma(a-b+1)
    lt = (n-1)*log(2*T+1) - log(8)
    lo, hi = 0, 1
    while lb(hi+n-2, n-1) <= lt: hi *= 2
    while hi-lo > 1:
        mid = (lo+hi)//2
        if lb(mid+n-2, n-1) <= lt: lo = mid
        else: hi = mid
    m = lo
    # exact correction
    while not ok(m): m -= 1
    while ok(m+1): m += 1
    return m

def w_of(n, r, T): return T*sqrt(8*(n-1)*log(4*r))

out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

pr("n r m 0.7n(n-1) m>=0.7n(n-1) 4w m>=4w ratio m/(2w) 0.12sqrt((n-1)/ln4r)")
for n in [55, 64, 100, 300, 1000, 2136, 2137, 2678, 2679, 4277, 4278]:
    m = exact_m(n, n)
    for (lab, r) in [("n", n), ("nlnn", int(n*log(n))+1), ("n2", n*n)]:
        w = w_of(n, r, n)
        pr(n, lab, r, m, round(0.7*n*(n-1),1), m >= 0.7*n*(n-1), round(4*w,1), m >= 4*w,
           round(m/(2*w),4), round(0.12*sqrt((n-1)/log(4*r)),4))

# sufficient condition n-1 >= 300 ln(4r) => m >= 4w : check at boundary r for several n
pr("--- boundary r = floor(exp((n-1)/300)/4) ---")
for n in [2100, 2500, 3000, 4000, 6000]:
    r = int(math.exp((n-1)/300)/4)
    if r < 1: continue
    m = exact_m(n, n); w = w_of(n, r, n)
    pr(n, r, m, round(4*w,1), m >= 4*w)

# Theorem 4.1 constant: W <= 8 ln|c|_1/m + 8 ln(e D C/(2 m)) with ln|c|_1 <= 3 ln U + 2(m-1) ln n,
# U <= (2n+1)^(n-1); claimed <= 13.46 + 16 ln n + 8 ln C  (< 14 + 16 ln n + 8 log+ C)
pr("--- Theorem 4.1 constant (C=1) ---")
for n in [64, 65, 100, 1000, 2137, 5000]:
    m = exact_m(n, n); D = 2*(n-1)*n
    lU = (n-1)*log(2*n+1)
    lc1 = 3*lU + 2*(m-1)*log(n)
    rhs = 8*lc1/m + 8*(log(D/2) - math.lgamma(m+1)/m)   # exact ln m! instead of m ln(m/e)
    rhs_paper = 8*lc1/m + 8*log(math.e*D/(2*m))
    pr(n, m, "rhs_exact_lnfact", round(rhs - 16*log(n),3), "rhs_with_mln(m/e)", round(rhs_paper-16*log(n),3),
       "claimed<=13.46", rhs_paper - 16*log(n) <= 13.46, "eD/2m", round(math.e*D/(2*m),4))

open(sys.argv[1] if len(sys.argv)>1 else "/dev/null","w").write("\n".join(out)+"\n")
