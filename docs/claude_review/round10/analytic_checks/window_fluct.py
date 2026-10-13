# Evidence for Remark 3.4(1): in the fake windows [P, P e^(1/(8q))] the genuine Gaussian-prime angle ensemble
# fluctuates (box counts, Hecke sums) by far more than the 2s = 4q points a graft moves.
import numpy as np, math
def run(P, q, L=64, K=2000):
    hi = int(P * math.exp(1.0/(8*q)))
    lo = P
    # segmented sieve on [lo, hi]
    seg = np.ones(hi - lo + 1, dtype=bool)
    for p in range(2, math.isqrt(hi) + 1):
        if p > 2 and p % 2 == 0: continue
        isp = all(p % d for d in range(2, math.isqrt(p) + 1))
        if not isp: continue
        start = max(p*p, ((lo + p - 1)//p)*p)
        seg[start - lo::p] = False
    angs = []
    for a in range(1, math.isqrt(hi) + 1):
        b_lo = math.isqrt(max(lo - a*a, 0)); b_hi = math.isqrt(max(hi - a*a, 0))
        if b_hi < 1: continue
        b = np.arange(max(b_lo, 1), b_hi + 1)
        n = a*a + b*b
        m = (n >= lo) & (n <= hi)
        b = b[m]; n = n[m]
        ok = seg[n - lo] & (n % 4 == 1)
        if ok.any(): angs.append(np.arctan2(b[ok], a))
    th = np.concatenate(angs)          # Gaussian primes in the open first quadrant, both conjugates
    n = len(th); s2 = 4*q              # graft moves 2s = 2*(2q) elements of G
    cnt = np.histogram(th, bins=L, range=(0, np.pi/2))[0]
    dev = np.abs(cnt - n/L)
    ks = np.arange(1, K + 1)
    S = np.abs(np.exp(4j*np.outer(ks, th)).sum(axis=1)) if n*K < 4e8 else None
    print("P=%.0e q=%d (M=%d): window [%d, %d], Gaussian primes n=%d, graft moves 2s=%d" % (P, q, q+1, lo, hi, n, s2))
    print("  %d equal sectors: mean count %.1f, max |dev| %.1f, median |dev| %.1f, sqrt(mean) %.1f"
          % (L, n/L, dev.max(), np.median(dev), math.sqrt(n/L)))
    if S is not None:
        print("  Hecke sums |S(k)|, 1<=k<=%d: median %.1f, min %.2f, max %.1f; sqrt(n)=%.1f; fraction of k with |S(k)| < 2s: %.3f"
              % (K, np.median(S), S.min(), S.max(), math.sqrt(n), np.mean(S < s2)))
for P, q in ((10**7, 11), (10**8, 11), (10**8, 23), (10**9, 47)):
    run(P, q)
