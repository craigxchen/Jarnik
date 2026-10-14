"""Strengthened J5 using projection-cheap Fano blocks (rho(P_j) <= 3j, proved in refJ.md §1.6).
Sylvester n = 2^k, C = 1.  Families: s Fano blocks {a,b,a+b} + c single columns, all on linearly independent bits:
image subset P_j^s x [-j,j]^c, rho <= (3s+2c) j, span = functions f_1(bits) + ... , 2^{2s+c} cells of size n/2^{2s+c}.
(A) uniform Q0: Q0 < Q4: L_4 point (Q0 abar >= 2 ln 15); Q4 <= Q0 < Q5: (s,c)=(2,1) point, cost 8j < m' (Q0 abar >= 2 ln 31);
    Q0 >= Q5: case B with L_1 (o_1 = m'-2j), needs o_1 Q5 > d, level g*(Q5) (g* decreasing).
(B) per-subspace radii: span (s,c), o = m' - (3s+2c) j, Q = N2/P, failure iff oQ > d and abar >= g*."""
import math
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 60
out = []
def P(s): print(s, flush=True); out.append(s)
def lncomb(a, b): return mp.loggamma(a + 1) - mp.loggamma(b + 1) - mp.loggamma(a - b + 1)
for k in range(12, 24):
    n = 2**k; d = n - 1
    w = mp.sqrt(8 * d * mp.log(4 * d)); j = int(mp.floor(w / 2)); m = int(mp.ceil(4 * w))
    assert min(abs(w/2 - mp.nint(w/2)), abs(4*w - mp.nint(4*w))) > mp.mpf(10)**-20
    if k <= 20: R1 = ("exact", 8 * math.comb(m + d - 1, d) <= 3**d)
    else: R1 = ("lgamma", mp.log(8) + lncomb(m + d - 1, d) <= d * mp.log(3))
    regime = d >= 316 * math.log(4 * d)
    N2 = Fr(n * (n - 1), 2)
    def Q(cells):
        s = n // cells; return N2 / (cells * Fr(s * (s - 1), 2))
    Q4, Q5 = Q(16), Q(32)
    assert 8 * j < m
    o1 = m - 2 * j; c = mp.mpf(n) / (n - 1)
    a = m * mp.log(4) - mp.log(mp.mpf(d) / (m + 1))
    q5 = mp.mpf(Q5.numerator) / Q5.denominator
    crit = o1 * q5 - d
    if crit > 0:
        gB = 2 * (a + o1 * mp.log(c * q5)) / crit
        sigma = mp.exp(-q5 * gB / 2) / 2; eta = 2 * sigma * q5 * c
        direct = m * mp.log(4) + o1 * mp.log(eta) <= mp.log(mp.mpf(d) / (m + 1)) - d * gB / 2 + mp.mpf(10)**-30
        side = eta <= 0.5 and (1 - eta) / (2 * (1 + eta)) >= sigma
    else:
        gB = mp.inf; direct = side = None
    gA1 = 2 * mp.log(15); gA2 = 2 * mp.log(31) / (mp.mpf(Q4.numerator) / Q4.denominator)
    gne = (3 * n * mp.log(n) - 14 * n) / (2 * d)
    P(f"k={k} n={n} j={j} m'={m} R1({R1[0]})={R1[1]} regime={regime} Q4={float(Q4):.5f} Q5={float(Q5):.5f} o1={o1} o1*Q5-d={mp.nstr(crit,8)}")
    P(f"   uniform-Q0 fails for every Q0>0 at every abar >= max(2ln15={mp.nstr(gA1,5)}, 2ln31/Q4={mp.nstr(gA2,4)}, g*(Q5)={mp.nstr(gB,6)}) ;"
      f" caseB direct={direct} side={side} ; near-extremal level {mp.nstr(gne,6)}")
    for (s_, c_) in [(1, 0), (2, 0), (1, 1), (1, 2), (0, 2), (0, 3)]:
        cells = 2**(2 * s_ + c_); o = m - (3 * s_ + 2 * c_) * j
        if o < 1: continue
        Qf = mp.mpf(Q(cells).numerator) / Q(cells).denominator; sz = n // cells
        Pf = cells * mp.mpf(sz * (sz - 1)) / 2
        cr = o * Qf - d
        if cr <= 0: continue
        K = mp.mpf(n) * sz / Pf
        gs = 2 * (m * mp.log(4) + o * mp.log(K) - mp.log(mp.mpf(d) / (m + 1))) / cr
        gmix = 2 * mp.log(cells - 1) + 1
        g = max(gs, gmix) + mp.mpf(10)**-40
        eta = K * mp.exp(-Qf * g / 2)
        direct = m * mp.log(4) + o * mp.log(eta) <= mp.log(mp.mpf(d) / (m + 1)) - d * g / 2
        mix = (mp.mpf(1) / (cells - 1) - eta) / (2 * (1 + eta)) >= mp.exp(-g / 2) / 2 and eta <= 0.5
        P(f"   per-subspace span (s,c)=({s_},{c_}) cells={cells} o={o} oQ-d={mp.nstr(cr,7)}: fails for abar >= {mp.nstr(g,7)}"
          f" (W/W_min = {mp.nstr(2*d*g/(3*n*mp.log(n)-14*n),4)}); direct={direct} mixing/eta_ok={mix}")
open(__file__.replace('.py', '_out.txt'), 'w').write("\n".join(out) + "\n")
