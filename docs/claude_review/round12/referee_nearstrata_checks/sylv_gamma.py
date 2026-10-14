# Gamma of equally spaced point u_t = t/(n-1), t=0..n-1 (binary-fraction point normalised):
# Gamma = (2/N2) sum_{delta=1}^{n-1} (n-delta) ln((n-1)/delta)
import mpmath as mp
mp.mp.dps = 40
def gam(n):
    N2 = mp.mpf(n)*(n-1)/2
    # sum (n-delta) ln(n-1) - sum (n-delta) ln delta
    S1 = N2*mp.log(n-1)
    # sum_{delta=1}^{n-1} (n-delta) ln delta = n ln((n-1)!) - sum delta ln delta
    # use exact-ish: loggamma and hyperfactorial: sum_{d=1}^{m} d ln d = ln H(m)
    m = n-1
    lnfact = mp.loggamma(m+1)
    # ln hyperfactorial via mpmath: hyperfac
    lnH = mp.log(mp.hyperfac(m)) if m < 5000 else None
    if lnH is None:
        # Euler-Maclaurin: sum_{d=1}^m d ln d = m^2/2 ln m - m^2/4 + m/2 ln m + 1/12 ln m + 1/12 - zeta'(-1)... use asymptotic of ln H(m)
        lnH = (mp.mpf(m)**2/2 + mp.mpf(m)/2 + mp.mpf(1)/12)*mp.log(m) - mp.mpf(m)**2/4 + mp.log(mp.glaisher) \
              + 1/(720*mp.mpf(m)**2) - 1/(5040*mp.mpf(m)**4) + 1/(10080*mp.mpf(m)**6)
    S2 = n*lnfact - lnH
    return 2*(S1-S2)/N2
for n in [8,64,1024,2**10+1, 2**12, 2**15, 2**16, 2**17, 2**20, 2**30, 2**60]:
    print(n, mp.nstr(gam(n), 15))
# direct check for small n
def gam_direct(n):
    N2 = mp.mpf(n)*(n-1)/2
    return 2*mp.fsum((n-dl)*mp.log(mp.mpf(n-1)/dl) for dl in range(1,n))/N2
for n in [8,64,1024,4096]:
    print('direct', n, mp.nstr(gam_direct(n),15))
