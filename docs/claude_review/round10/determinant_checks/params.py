# Parameter feasibility for Theorem 3 (slab-box) : exact integer binomials, logs via mpmath-free math.
import math
def lbinom(a,b): return math.lgamma(a+1)-math.lgamma(b+1)-math.lgamma(a-b+1)
def params(n,r,T):
    w = T*math.sqrt(8*(n-1)*math.log(4*r))
    lU0 = (n-1)*math.log(2*T+1)          # ln |box|
    # m maximal with binom(m+n-2,n-1) <= (2T+1)^{n-1}/8   (so U >= |box|/2 >= 4E)
    lo,hi=0,1
    while lbinom(hi+n-2,n-1) <= lU0-math.log(8): hi*=2
    while hi-lo>1:
        mid=(lo+hi)//2
        if lbinom(mid+n-2,n-1) <= lU0-math.log(8): lo=mid
        else: hi=mid
    return w, lo
def W1(n,T,m,w,C):
    # W threshold: (m/2 - w) W > 2 ln|c|_1 + 2 m ln(D C/2) - 2 ln m!
    D = 2*(n-1)*T
    lU = (n-1)*math.log(2*T+1); lB = (m-1)*math.log(T)
    lH = 2*(lU+lB); lc1 = lU + lH
    rhs = 2*lc1 + 2*m*math.log(D*C/2) - 2*math.lgamma(m+1)
    return rhs/(m/2-w) if m/2 > w else float('inf')
print("n  r   T   m   4w   ok  W1(C=1)  W1(C=10)  ratio m/(2w)")
for n in [64,128,256,512,1024,2048,4096]:
    for r in [n, int(n*math.log(n))]:
        T = n
        w,m = params(n,r,T)
        print(n, r, T, m, round(4*w,1), m>=4*w, round(W1(n,T,m,w,1.0),2), round(W1(n,T,m,w,10.0),2), round(m/(2*w),3))
# smallest n (T=n) with m >= 4w for r = n and r = n ln n
for label,f in [("r=n",lambda n:n),("r=n ln n",lambda n:int(n*math.log(n))+1),("r=n^2",lambda n:n*n)]:
    lo,hi=16,1<<16
    while hi-lo>1:
        mid=(lo+hi)//2
        w,m=params(mid,f(mid),mid)
        if m>=4*w: hi=mid
        else: lo=mid
    n=hi
    print("threshold", label, "n0 =", n, " 300 ln(4r) =", round(300*math.log(4*f(n)),1))
# check the closed-form lower bound m >= (n-1)(2T+1)/(e 8^{1/(n-1)}) - (n-1)
for n in [20,50,100,400]:
    for T in [1,5,n]:
        w,m = params(n,n,T)
        lb = (n-1)*(2*T+1)/(math.e*8**(1/(n-1))) - (n-1)
        assert m >= lb, (n,T,m,lb)
print("closed-form lower bound for m verified on grid")
