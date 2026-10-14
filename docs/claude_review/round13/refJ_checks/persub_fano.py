"""Per-subspace collision hypothesis ZE-TC^coll_sigma on Sylvester boxes, using the 3-dim Fano span
U = span(1_{S_a}, 1_{S_b}, 1_{S_{a+b}}) (4 cells = cosets of {a.x = b.x = 0}).  Proved: pi_U(A) subset P_j, rho(P_j) <= 3j,
so transversal order o_U >= m' - 3j (Lemma 5.0).  Same test point and distance argument as for the coordinate span L_2
(same 4-cell collision partition, P(U) = P(L_2)).  Failure iff o Q_2 > d and abar >= g* (formula as in j5_recheck).
Also: conditional line with o = m' - 4j for the 7-dim XOR family of a 3-dim V (8 cells), rho <= 4j being only EVIDENCE (j<=2)."""
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 60
out = []
def P(s): print(s); out.append(s)
for k in range(12, 19):
    n = 2**k; d = n - 1
    w = mp.sqrt(8 * d * mp.log(4 * d)); j = int(mp.floor(w / 2)); m = int(mp.ceil(4 * w))
    N2 = Fr(n * (n - 1), 2)
    for label, b, o in [("Fano span U (proved rho<=3j)", 2, m - 3 * j), ("coordinate L_2 (rho<=4j)", 2, m - 4 * j),
                        ("8-cell XOR span (CONDITIONAL rho<=4j)", 3, m - 4 * j), ("coordinate L_3 (rho<=6j)", 3, m - 6 * j)]:
        s = n >> b; Pb = 2**b * Fr(s * (s - 1), 2); Qb = N2 / Pb
        Qf = mp.mpf(Qb.numerator) / Qb.denominator; Pf = mp.mpf(Pb.numerator) / Pb.denominator
        crit = o * Qf - d
        if crit <= 0:
            P(f"k={k} {label}: o={o} oQ-d={mp.nstr(crit,7)} <= 0: no failure"); continue
        K = mp.mpf(n) * s / Pf
        gs = 2 * (m * mp.log(4) + o * mp.log(K) - mp.log(mp.mpf(d) / (m + 1))) / crit
        g = gs + mp.mpf(10)**-40
        eta = K * mp.exp(-Qf * g / 2)
        direct = (m * mp.log(4) + o * mp.log(eta)) <= (mp.log(mp.mpf(d) / (m + 1)) - d * g / 2)
        mix = (mp.mpf(1) / (2**b - 1) - eta) / (2 * (1 + eta)) >= mp.exp(-g / 2) / 2
        Wmin = 3 * n * mp.log(n) - 14 * n
        P(f"k={k} {label}: o={o} oQ-d={mp.nstr(crit,7)} fails for abar >= {mp.nstr(gs,8)} (W >= {mp.nstr(2*d*gs,8)} = {mp.nstr(2*d*gs/Wmin,4)} W_min); direct={direct} eta={mp.nstr(eta,3)} mixing_ok={mix}")
open(__file__.replace('.py', '_out.txt'), 'w').write("\n".join(out) + "\n")
