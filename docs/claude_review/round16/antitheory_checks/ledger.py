"""Loss ledger (C = 1, W = W_min, m = ceil(4w)). All logs natural.
L0 = ln(1/eps_0)                      : Step-0 budget
Lc = m ln(2 d e/m)                    : base cost of an exponent-m certificate (m!/(2d)^m)
Lpoly = m ln n                        : extra cost at distance 1/n from the strata, exponent m
L16 = m (8 abar + ln 2)               : extra cost at the separation radius of a 16-class cheap span, exponent m
Edg = needed Edgeworth/LLT precision exponent: precision e^{-(Lc + Lpoly)}"""
from common import *
print("n r m L0 Lc Lc+Lpoly room_poly=L0-Lc-Lpoly Lc+L16 room_16 (Lc+Lpoly)/ln(n)")
for n, r in [(2900, 2900), (4096, 4095), (10**4, 10**4), (65536, 65535), (10**5, 10**5), (473716, 473716), (10**6, 10**6), (10**9, 10**9), (10**12, 10**12)]:
    d, w, j, m = wjm(n, r)
    W = Wmin(n); abar = W / (2 * d)
    L0 = ln_inv_eps0(n, m, W)
    Lc = m * log(2 * d * mp.e / m)
    Lp = m * log(n)
    L16 = m * (8 * abar + log(2))
    print(n, r, m, mp.nstr(L0, 7), mp.nstr(Lc, 7), mp.nstr(Lc + Lp, 7), mp.nstr(L0 - Lc - Lp, 7), mp.nstr(Lc + L16, 7), mp.nstr(L0 - Lc - L16, 7), mp.nstr((Lc + Lp) / log(n), 6))
