"""J5b: per-subspace collision hypothesis (radius sigma_pi = (1/2) exp(-N2 g/(2 P(pi))) for EVERY collision subspace,
as Lemma 6.6 / Theorem 6.8 (first part) actually supply) on Sylvester boxes, coordinate spans L_b (b columns).
Test point u = y + eta v, y = class value c/(2^b-1) (base class 0), v = rank/s within each class (s = n/2^b, base rank 0).
Refinements pi of pi_b: dist(u, Delta_pi) >= eta (2P(pi)/n)/(2s); need >= 2 sigma_pi (normalisation |u|_inf <= 2);
f(P) = (n s/(2P)) e^{-N2 g/(2P)} is increasing in P for P < N2 g/2, so the binding pi is pi_b:
eta = (n s/P_b) e^{-Q_b g/2}. Partitions mixing classes: dist >= (1/(2^b-1) - 2 eta)/(2(1+eta)) >= sigma_pi as soon as
g >= 2 ln(2^b - 1) + 1 (checked below). Violation iff 4^{m'} eta^{o_b} <= (d/(m'+1)) e^{-d g/2}:
g >= gstar_b = 2[m' ln 4 + o_b ln(n s/P_b) - ln(d/(m'+1))]/(o_b Q_b - d)  when o_b Q_b > d.
"""
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 50
out = []
def P(s): print(s); out.append(s)
for k in range(12, 15):
    n = 2**k; d = n - 1
    w = mp.sqrt(8 * d * mp.log(4 * d)); j = int(mp.floor(w / 2)); m = int(mp.ceil(4 * w))
    bmax = (m - 1) // (2 * j); N2 = Fr(n * (n - 1), 2)
    for b in range(1, bmax + 1):
        s = n // 2**b; Pb = 2**b * Fr(s * (s - 1), 2); Qb = N2 / Pb; o = m - 2 * j * b
        Qf = mp.mpf(Qb.numerator) / Qb.denominator
        crit = o * Qf - d
        if crit > 0:
            gs = 2 * (m * mp.log(4) + o * mp.log(mp.mpf(n) * s / (mp.mpf(Pb.numerator) / Pb.denominator)) - mp.log(mp.mpf(d) / (m + 1))) / crit
            g = max(gs, 2 * mp.log(2**b - 1) + 1)
            eta = (mp.mpf(n) * s / (mp.mpf(Pb.numerator) / Pb.denominator)) * mp.exp(-Qf * g / 2)
            mixok = (mp.mpf(1) / (2**b - 1) - 2 * eta) / (2 * (1 + eta)) >= mp.mpf(1) / 2 * mp.exp(-g / 2)
            P(f"k={k} b={b} o_b={o} Q_b={float(Qb):.5f} o_bQ_b-d={mp.nstr(crit,6)} gstar={mp.nstr(gs,6)} fail_from_g={mp.nstr(g,6)} eta={mp.nstr(eta,4)} mixed_ok={mixok} W_fail={mp.nstr(2*d*g,6)} W_min={mp.nstr(3*n*mp.log(n)-mp.mpf('13.91')*n,6)}")
        else:
            P(f"k={k} b={b} o_b={o} Q_b={float(Qb):.5f} o_bQ_b-d={mp.nstr(crit,6)} : no violation from this span (proved order)")
open(__file__.replace('.py', '_out.txt'), 'w').write("\n".join(out) + "\n")
