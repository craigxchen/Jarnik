"""Referee recomputation of J5 (round12/README.md §4.1), independent code.

Sylvester profile, n = 2^k, d = n-1, r = d sign columns, C = 1 unless stated (C' = 1), abar = g = W/(2d),
eps0 = C (d/(m'+1)) e^{-d g/2}.  m' = ceil(4w), j = floor(w/2), w = sqrt(8 d ln 4r).

For each k we
 (1) certify floor/ceil (distance of w/2, 4w from integers at 80 digits),
 (2) check R1: 8 binom(m'+d-1, d) <= 3^d with exact integers (k <= 17; k>=18 via lgamma, evidence),
 (3) compute P_b, Q_b exactly (Fractions), o_b = m' - 2jb (b <= 4; o_4 >= 1 checked),
 (4) case A threshold Q0*g >= 2 ln 15 (worst Q0 = 1),
 (5) case B: the failure inequality  m' ln4 + o1 ln eta <= ln C + ln(d/(m'+1)) - d g/2  with
     eta = 2 sigma Q0 n/(n-1), sigma = e^{-Q0 g/2}/2, checked DIRECTLY at g = g*(Q4) (+1e-30) and on a
     dense Q0 grid in [Q4, 1e6] plus the analytic monotonicity; side conditions eta <= 1/2 and the
     mixing distance (1-eta)/(2(1+eta)) >= sigma also checked directly,
 (6) C-dependence: extra term 2 ln(1/C)/(o1 Q4 - d),
 (7) per-subspace radii, coordinate spans b = 1..4 (own derivation of the binding refinement).
"""
import math
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 80
out = []
def P(s=""):
    print(s); out.append(s)

def lncomb(a, b):
    return mp.loggamma(a + 1) - mp.loggamma(b + 1) - mp.loggamma(a - b + 1)

def frac(x):
    return mp.mpf(x.numerator) / x.denominator

for k in range(12, 21):
    n = 2**k; d = n - 1; r = d
    w = mp.sqrt(8 * d * mp.log(4 * r))
    j = int(mp.floor(w / 2)); m = int(mp.ceil(4 * w))
    dist_int = min(abs(w / 2 - mp.nint(w / 2)), abs(4 * w - mp.nint(4 * w)))
    assert dist_int > mp.mpf(10)**-20
    if k <= 17:
        R1 = 8 * math.comb(m + d - 1, d) <= 3**d
        R1s = f"R1(exact)={R1}"
    else:
        R1 = mp.log(8) + lncomb(m + d - 1, d) <= d * mp.log(3)
        R1s = f"R1(lgamma,evidence)={R1}"
    # m_max by lgamma bisection (evidence only; only m' <= m_max is used, and that is exact above)
    lo, hi = 1, d
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mp.log(8) + lncomb(mid + d - 1, d) <= d * mp.log(3): lo = mid
        else: hi = mid - 1
    mmax = lo
    N2 = Fr(n * (n - 1), 2)
    Pb = {b: 2**b * Fr((n >> b) * ((n >> b) - 1), 2) for b in range(1, 6)}
    Qb = {b: N2 / Pb[b] for b in Pb}
    ob = {b: m - 2 * j * b for b in Pb}
    assert ob[4] >= 1 and ob[5] <= 0
    regime = (d >= 316 * math.log(4 * r))
    Q4 = frac(Qb[4]); o1 = ob[1]
    c = mp.mpf(n) / (n - 1)
    a = m * mp.log(4) - mp.log(mp.mpf(d) / (m + 1))
    def gstar(Q, C=1):
        den = o1 * Q - d
        if den <= 0: return mp.inf
        return 2 * (a + o1 * mp.log(c * Q) + mp.log(1 / mp.mpf(C))) / den
    def direct_fail(Q, g, C=1):
        sigma = mp.exp(-Q * g / 2) / 2
        eta = 2 * sigma * Q * c
        lhs = m * mp.log(4) + o1 * mp.log(eta)
        rhs = mp.log(C) + mp.log(mp.mpf(d) / (m + 1)) - d * g / 2
        side = (eta <= mp.mpf(1) / 2) and ((1 - eta) / (2 * (1 + eta)) >= sigma) and (mp.mpf(1) / 30 >= 0)
        return lhs <= rhs, side, eta
    gB = gstar(Q4)
    if gB < mp.inf:
        ok, side, eta = direct_fail(Q4, gB + mp.mpf(10)**-40)
        # dense grid: failure at g = max over grid of gstar must hold for each Q
        grid_ok = True; worst = mp.mpf(0)
        Q = Q4
        while Q < 10**6:
            gq = gstar(Q)
            if gq > gB + mp.mpf(10)**-40: grid_ok = False
            f, s, _ = direct_fail(Q, gB + mp.mpf(10)**-40)
            if not (f and s): grid_ok = False
            Q *= mp.mpf('1.01')
        Cdep = 2 / (o1 * Q4 - d)
    else:
        ok = side = grid_ok = None; eta = None; Cdep = None
    gA = 2 * mp.log(15)
    gne = (3 * n * mp.log(n) - 14 * n) / (2 * d)
    P(f"k={k} n={n} w={mp.nstr(w,10)} j={j} m'={m} m_max~{mmax} {R1s} regime(d>=316ln4r)={regime}")
    P(f"   Q_b={[mp.nstr(frac(Qb[b]),8) for b in range(1,5)]} o_b={[ob[b] for b in range(1,5)]}")
    P(f"   caseA: fails for Q0<Q4 when Q0*g >= 2ln15={mp.nstr(gA,6)};  caseB: o1*Q4-d={mp.nstr(o1*Q4-d,8)} g*(Q4)={mp.nstr(gB,8)}"
      f" direct_check={ok} side_conditions={side} eta={mp.nstr(eta,5) if eta is not None else '-'} grid[Q4,1e6]_ok={grid_ok}")
    P(f"   C-dependence of case B threshold: + {mp.nstr(Cdep,5) if Cdep else '-'} * ln(1/C);  near-extremal level g_ne=(3n ln n-14n)/(2d)={mp.nstr(gne,6)};"
      f" uniform failure for all Q0 at g >= {mp.nstr(max(gA, gB),6)}")
    # per-subspace: refinements pi of pi_b with P(pi) = P: block of size >= 1 + 2P/n inside a class (size s);
    # v = rank/s in each class: dist(u,Delta_pi) >= eta*(2P/n)/(2s) = eta P/(n s). Need >= (1+eta) sigma_pi.
    # sigma_pi/(P) with sigma_pi = e^{-N2 g/(2P)}/2 : need eta >= (1+eta) n s e^{-N2 g/(2P)}/(2P) for all 1<=P<=P_b.
    # f(P) = e^{-N2 g/(2P)}/P increasing on P < N2 g/2, and P_b <= N2 < N2 g/2 for g > 2, so binding at P = P_b.
    for b in range(1, 5):
        s = n >> b; Pbf = frac(Pb[b]); Qf = frac(Qb[b]); o = ob[b]
        crit = o * Qf - d
        if crit <= 0:
            P(f"   per-subspace b={b}: o_b Q_b - d = {mp.nstr(crit,8)} <= 0 : no failure from this span"); continue
        # eta = 2 n s/(2 P_b) e^{-Q_b g/2} * (1+eta) ~ ; take eta = 2*(n s/(2 Pb)) e^{-Qb g/2} (factor 2 covers 1+eta<=2)
        K = mp.mpf(n) * s / Pbf
        gs = 2 * (m * mp.log(4) + o * mp.log(K) - mp.log(mp.mpf(d) / (m + 1))) / crit
        g = gs + mp.mpf(10)**-40
        eta = K * mp.exp(-Qf * g / 2)
        lhs = m * mp.log(4) + o * mp.log(eta); rhs = mp.log(mp.mpf(d) / (m + 1)) - d * g / 2
        # mixing partitions: values y on 2^b classes are (N(x)-N(base))/(2^b-1) with gaps >= 1/(2^b-1),
        # a mixing block has spread >= 1/(2^b-1) - eta, dist >= that/2, normalised by <= 1+eta; sigma_pi <= e^{-g/2}/2.
        mix = (mp.mpf(1) / (2**b - 1) - eta) / (2 * (1 + eta)) >= mp.exp(-g / 2) / 2
        P(f"   per-subspace b={b}: o_b Q_b - d = {mp.nstr(crit,8)} > 0, fails for g >= {mp.nstr(gs,8)} (W >= {mp.nstr(2*d*gs,7)},"
          f" W/W_min = {mp.nstr(2*d*gs/(3*n*mp.log(n)-14*n),4)}); direct={lhs<=rhs} eta={mp.nstr(eta,4)}<=1/2:{eta<=0.5} mixing_ok={mix}")
open(__file__.replace('.py', '_out.txt'), 'w').write("\n".join(out) + "\n")
