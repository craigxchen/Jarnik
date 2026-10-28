"""Independent recomputation for referee_ze.md.
Exact big-integer binomials; logs via math.lgamma only for display.
(a) m_max, w, ceil(4w), W_min, E(m,1), J_lemma (Cor 3.2).
(b) Prop 3.5(1) margin, GV S_d, Prop 4.4 margins.
(c) Thm 3.6 barrier as a function of lambda (threshold lambda*), mu = ceil(4w), m_max.
(d) Thm 3.6 with the AGL(k,2) automorphism orbit of the Sylvester profile (n = 2^k)."""
import math
from fractions import Fraction

def lnb(a, b):
    return math.lgamma(a+1)-math.lgamma(b+1)-math.lgamma(a-b+1)

def m_max(d):
    t = 3**d; lo, hi = 0, 2*d
    while lo < hi:
        mid = (lo+hi+1)//2
        if 8*math.comb(mid+d-1, d) <= t: lo = mid
        else: hi = mid-1
    return lo

def tau(mu, j, lam, n, W):
    # floor(mu(1 + 4ln2/W - j lam/sqrt n)) + 1 ; float is fine away from integer boundaries (checked below)
    x = mu*(1+4*math.log(2)/W - j*lam/math.sqrt(n))
    if abs(x-round(x)) < 1e-6: print("WARNING near-integer", mu, j, x)
    return math.floor(x)+1

def count(mu, lam, n, W):
    d = n-1; J = math.floor(1.07*math.sqrt(n)/lam); c = 0
    for j in range(1, J+1):
        t = tau(mu, j, lam, n, W)
        if t <= 0: break
        c += 2*math.comb(t+d-1, d)
    return c

def avail(mu, d):
    return (3**d+1)//2 - math.comb(mu+d-1, d)

def ln(x):
    return math.log(x) if x > 0 else float('-inf')

def agl_order(k):
    g = 1
    for i in range(k): g *= (2**k - 2**i)
    return (2**k)*g

out = []
for n in [2966, 3000, 4096, 8192]:
    d = n-1; W = 3*n*math.log(n)-14*n
    w = math.sqrt(8*d*math.log(4*n)); m4 = math.ceil(4*w); mm = m_max(d)
    # (a) E and J_lemma
    def E(m, j):
        lnc = 3*d*math.log(3)
        return (2*lnc + 2*m*math.log(2*d*j/2) - 2*math.lgamma(m+1))/(W/2)
    def Jl(m):
        j = 0
        while m-2*(j+1)*w > E(m, j+1): j += 1
        return j
    out.append(f"n={n} d={d} W={W:.1f} w={w:.2f} ceil4w={m4} m_max={mm} E(mmax,1)={E(mm,1):.3f} Jl(mmax)={Jl(mm)} Jl(4w)={Jl(m4)}")
    # (b)
    lnd2 = math.lgamma(d+1)+d*math.log(2)
    out.append(f"   Prop3.5(1): ln(d!2^d)={lnd2:.1f}; minus ln binom(1+d,d)={lnd2-math.log(d+1):.1f}")
    Sd = lnd2 - lnb(3*d, d)
    out.append(f"   GV S_d={Sd:.1f}; log10 W bound={(Sd-math.log(2.14*d/math.log(5)))/math.log(10):.1f}")
    l44 = math.log(d+1)+d*math.log(m4-1)
    out.append(f"   Prop4.4: ln((d+1)(m-1)^d)={l44:.1f}, minus d ln3 = {l44-d*math.log(3):.1f}")
    # (c) lambda threshold
    lamS = n/(n-1)
    for mu in (m4, mm):
        c = count(mu, lamS, n, W); a = avail(mu, d)
        # threshold lambda via log-sum-exp (float; evidence), then exact check either side
        la = ln(a)
        def lcount(lam):
            J = math.floor(1.07*math.sqrt(n)/lam); terms = []
            for j in range(1, J+1):
                t = math.floor(mu*(1+4*math.log(2)/W - j*lam/math.sqrt(n)))+1
                if t <= 0: break
                terms.append(math.log(2)+lnb(t+d-1, d))
            mx = max(terms); return mx+math.log(sum(math.exp(x-mx) for x in terms))
        lo, hi = 1e-3, lamS
        if lcount(lo) < la: thr = "<=1e-3"
        else:
            for _ in range(50):
                mid = (lo+hi)/2
                if lcount(mid) < la: hi = mid
                else: lo = mid
            ex_hi = count(mu, hi*1.01, n, W) < a; ex_lo = count(mu, lo*0.99, n, W) < a
            thr = f"{hi:.5f} [exact: holds at 1.01*thr: {ex_hi}; at 0.99*thr: {ex_lo}]"
        out.append(f"   mu={mu}: lam=n/(n-1): ln count={ln(c):.1f} ln avail={ln(a):.1f} holds={c<a}; barrier holds iff lam >~ {thr} (2.14 sqrt n/d = {2.14*math.sqrt(n)/d:.4f})")
    # (d) Sylvester automorphism orbit
    k = round(math.log2(n))
    if 2**k == n:
        g = agl_order(k)
        for mu in (m4, mm):
            c = count(mu, lamS, n, W)*g; a = avail(mu, d)
            out.append(f"   AGL({k},2) orbit: ln|AGL|={math.log(g):.2f}; mu={mu}: ln(count*|AGL|)={ln(c):.1f} ln avail={ln(a):.1f} holds={c<a}")
        # largest mu for which the orbit-inclusive count holds
        mu = m4
        while count(mu+1, lamS, n, W)*g < avail(mu+1, d): mu += 1
        out.append(f"   AGL orbit: largest mu with barrier = {mu} (m_max={mm})")
    print("\n".join(out), flush=True); out.clear()
