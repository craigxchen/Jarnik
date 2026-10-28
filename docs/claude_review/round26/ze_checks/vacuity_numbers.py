"""Numbers for Sections 3-5 of ze.md (T = 1 slab box, regime sizes).
Exact: m_max (big-integer binomials), comparisons of big integers where stated.
Floating (labelled): E(c,j) budget at W_min, logs of huge quantities (mpmath, 50 digits)."""
import math
from mpmath import mp, mpf, log, loggamma, sqrt, ceil, floor, exp
mp.dps=50
def lbinom(a,b): return loggamma(a+1)-loggamma(b+1)-loggamma(a-b+1)
def m_max_exact(d):
    # largest m with 8*binom(m+d-1,d) <= 3^d (exact integers)
    lo,hi=0,2*d
    t=3**d
    while lo<hi:
        mid=(lo+hi+1)//2
        if 8*math.comb(mid+d-1,d)<=t: lo=mid
        else: hi=mid-1
    return lo
print("T=1 slab box; r = n (and Sylvester r = n-1 for lam); C = 1; W = W_min = 3 n ln n - 14 n")
print(f"{'n':>7} {'m_max':>6} {'w':>8} {'4w':>6} {'E(m_max,1)':>10} {'J_lemma(m_max)':>14} {'J_lemma(ceil4w)':>15} {'sqrt(n)/lam_Syl':>15} {'d/2':>6}")
rows=[]
for n in [2966,3000,4096,8192,16384,10**5,10**6]:
    d=n-1; r=n
    w=sqrt(8*d*log(4*r))
    m=m_max_exact(d) if n<=16384 else None
    if m is None:
        # use proved bound m_max >= 0.6417 d - 2.2137 (t1rigour) as a floor (label: proved lower bound)
        m=int(floor(mpf('0.6417')*d-mpf('2.2137')))
    W=3*n*log(n)-14*n
    lnU=d*log(3)                      # U <= 3^d
    def E(mm,j):
        # 4 ln|c|_1 + 4 mm ln(D j C/2) - 4 ln mm!  over W,  |c|_1 <= U^3, D = 2d
        return (4*3*lnU + 4*mm*log(d*j) - 4*loggamma(mm+1))/W
    def J_lemma(mm):
        # largest j with mm - j*2w > E(mm,j)  (Lam(c) <= h_A-bound 2w for every c on A)
        j=0
        while mm-(j+1)*2*w>E(mm,j+1): j+=1
        return j
    m4=int(ceil(4*w))
    lamS=mpf(n)/(n-1)
    rows.append((n,m,float(w),m4))
    print(f"{n:>7} {m:>6} {float(w):8.2f} {m4:>6} {float(E(m,1)):10.3f} {J_lemma(m):>14} {J_lemma(m4):>15} {float(sqrt(n)/lamS):15.2f} {d//2:>6}")
print()
print("Philippon (torus picture), degree normalisation d!(2T)^d at T=1; singleton Sigma, T0 <= mu/d <= 2:")
for n,m,w,m4 in rows:
    d=n-1
    lhs=lbinom(2+d,d)            # binom(T0+d,d) with T0<=2
    rhs=loggamma(d+1)+d*log(2)
    print(f"  n={n}: ln binom(T0+d,d) <= {float(lhs):.1f}  vs  ln d!2^d = {float(rhs):.1f}  -> H={{e}} inequality holds (vacuous)")
print()
print("Generous variant (vanishing only on Sigma, order T0 <= mu <= 2d, |Sigma| <= 2J*+1, J* <= T d W/log 5):")
for n,m,w,m4 in rows:
    d=n-1
    lhs_noJ=lbinom(3*d,d)
    rhs=loggamma(d+1)+d*log(2)
    # vacuous iff ln(2J*+1) <= rhs - lhs_noJ ; J* <= d W/log5  => vacuous for W <= W_vac
    slack=rhs-lhs_noJ
    print(f"  n={n}: ln[d!2^d / binom(3d,d)] = {float(slack):.1f}; vacuous for all W with ln(2dW/log5+1) <= that, i.e. log10 W_vac ~ {float((slack-log(2*d/log(5)))/log(10)):.0f}")
print()
print("Dual picture G_a^d x G_m, bidegree (m-1,1): |Sigma| <= |A| <= 3^d  vs  (d+1)(m-1)^d")
for n,m,w,m4 in rows:
    d=n-1
    for mm in (m4,m):
        print(f"  n={n} m={mm}: ln 3^d = {float(d*log(3)):.1f} <= ln (d+1)(m-1)^d = {float(log(d+1)+d*log(mm-1)):.1f}  (m-1 >= 3 suffices)")
print()
print("Single-polynomial count (Thm 4.1): conditions sum_{1<=|j|<=J} binom(tau_j+d-1,d), tau_j = ceil(mu(1-|j|lam/sqrt n)), vs dim V_H >= 3^d/4")
for n,m,w,m4 in rows[:5]:
    d=n-1; lam=mpf(n)/(n-1)
    for mu in (m4,m):
        tot=mpf(0); J=int(floor(sqrt(n)/lam))
        for j in range(1,J+1):
            tau=int(ceil(mu*(1-j*lam/sqrt(n))))
            if tau<=0: break
            tot+=2*exp(lbinom(tau+d-1,d))
        print(f"  n={n} mu={mu}: ln(conditions) = {float(log(tot)) if tot>0 else float('-inf'):.1f}   ln(3^d/4) = {float(d*log(3)-log(4)):.1f}")
