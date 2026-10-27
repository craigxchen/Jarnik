# Skeptic check of steelman S3: Gaussian-proxy Carbery-Wright bound for AC(u).
# Pr(|Q|<=delta) <= C_CW * m' * (delta/sigma)^(1/m'),  sigma >= sqrt(m'!) (|u|_2>=1),
# delta = d^m' * eps_0,  eps_0 = C0 * d * exp(-W/4)/(m'+1)   (record round12 notation).
# AC needs the bound < 1/2.  Report log10 of the bound for several (n, r, W).
import math
def logbound(n, r, W, C_CW=1.0, C0=1.0):
    d = n - 1
    w = math.sqrt(8*d*math.log(4*r))
    m = math.ceil(4*w)
    ln_eps0 = math.log(C0*d/(m+1)) - W/4
    ln_delta = m*math.log(d) + ln_eps0
    ln_sigma = 0.5*math.lgamma(m+1)
    ln_pr = math.log(C_CW*m) + (ln_delta - ln_sigma)/m
    return m, W/(4*m), ln_pr/math.log(10)
print("n, r-label, W-label, m', W/(4m'), log10(CW bound)  [AC needs < log10(0.5) = -0.30]")
for n in [2900, 5809, 10**4, 3*10**4, 10**5, 10**6]:
    for rl, r in [("r=n", n), ("r=n^2", n*n)]:
        for wl, W in [("3nlnn", 3*n*math.log(n)), ("10nlnn", 10*n*math.log(n))]:
            m, ratio, lp = logbound(n, r, W)
            print(f"{n:8d} {rl:6s} {wl:7s} m'={m:6d} W/4m'={ratio:8.2f} log10Pr<={lp:9.2f}")
# crossover at W = 3 n ln n, r = n
for r_of in [("r=n", lambda n: n), ("r=n^2", lambda n: n*n)]:
    n = 100
    while logbound(n, r_of[1](n), 3*n*math.log(n))[2] > math.log10(0.5):
        n = int(n*1.05)+1
    print("crossover (C_CW=1, W=3 n ln n,", r_of[0], "): n ~", n)
