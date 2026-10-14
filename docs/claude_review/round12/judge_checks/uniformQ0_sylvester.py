"""J5: the uniform-Q0 collision hypothesis ZE-TC^coll_eps(A, ceil(4w); Q0) fails on Sylvester regime boxes.

Setting (t1rigour.md / nearstrata.md): T = 1, n = 2^k points x in F_2^k, Sylvester columns 1[a.x=0] (a != 0),
r = d = n-1, w = sqrt(8 d ln 4r), j = floor(w/2), m' = ceil(4w), C = C' = 1, abar = g = W/(2d),
eps0 = d e^{-W/4}/(m'+1) = (d/(m'+1)) e^{-d g/2}.
Coordinate columns S_i = {x : x_i = 0}: on A, |k(S_i)| <= j (sign-column slab), so b of them cost 2jb and the
image of A under them has <= (2j+1)^b points, interpolation degree <= 2jb; transversal order o_b >= m' - 2jb
(nearstrata Lemma 5.0: kappa(hat(y + eta v)) <= 4^{m'} eta^{o_b}, |y|_inf=1, |v|_inf<=1, eta<=1/2).
Exact: R1 (8 binom(m'+d-1,d) <= 3^d), cheapness 2jb <= m'-1, P_b, Q_b = N2/P_b (Fractions).
Failure levels: case A (Q0 < Q_bmax): g >= 2 ln(2^bmax - 1) (worst Q0 = 1); case B (Q0 >= Q_bmax):
g >= gstar(Q_bmax), gstar(Q) = 2[m' ln4 + o1 ln(Q n/(n-1)) - ln(d/(m'+1))]/(o1 Q - d), decreasing for Q >= e.
Also the per-subspace version (radius (1/2)exp(-N2 g/(2P(L)))) for each coordinate span: failure iff o_b Q_b > d.
"""
import math
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 60

def params(k):
    n = 2**k; d = n - 1; r = d
    w = mp.sqrt(8 * d * mp.log(4 * r))
    j = int(mp.floor(w / 2)); mprime = int(mp.ceil(4 * w))
    return n, d, r, w, j, mprime

def R1_exact(d, mprime):
    return 8 * math.comb(mprime + d - 1, d) <= 3**d

out = []
def P(s):
    print(s); out.append(s)

P("k n w j m' R1 bmax | Q_bmax | caseA_g caseB_gstar eta_at_gstar | g_nearextremal W1 | per-subspace o_b*Q_b - d (b=1..bmax)")
for k in range(12, 21):
    n, d, r, w, j, mp_ = params(k)
    r1 = R1_exact(d, mp_) if k <= 17 else None
    bmax = (mp_ - 1) // (2 * j)
    N2 = Fr(n * (n - 1), 2)
    Pb = {b: (2**b) * Fr((n // 2**b) * (n // 2**b - 1), 2) for b in range(1, bmax + 1)}
    Qb = {b: N2 / Pb[b] for b in Pb}
    ob = {b: mp_ - 2 * j * b for b in Pb}
    o1 = ob[1]
    Qm = Qb[bmax]
    gA = 2 * mp.log(2**bmax - 1)
    if o1 * Qm > d:
        num = mp_ * mp.log(4) + o1 * mp.log(mp.mpf(Qm.numerator) / Qm.denominator * n / (n - 1)) - mp.log(mp.mpf(d) / (mp_ + 1))
        gstar = 2 * num / (o1 * (mp.mpf(Qm.numerator) / Qm.denominator) - d)
        sigma = mp.mpf(1) / 2 * mp.exp(-(mp.mpf(Qm.numerator) / Qm.denominator) * gstar / 2)
        eta = 2 * sigma * (mp.mpf(Qm.numerator) / Qm.denominator) * n / (n - 1)
    else:
        gstar = mp.inf; eta = None
    gne = (3 * n * mp.log(n) - mp.mpf('13.91') * n) / (2 * d)
    W1 = 24 * d * mp.log(3) / mp_ + 8 * mp.log(mp.e * d / mp_)
    per = [int(ob[b] * Qb[b] - d) if True else None for b in range(1, bmax + 1)]
    P(f"{k} {n} {mp.nstr(w,8)} {j} {mp_} {r1} {bmax} | {float(Qm):.5f} | {mp.nstr(gA,5)} {mp.nstr(gstar,6)} {mp.nstr(eta,4) if eta is not None else '-'} | {mp.nstr(gne,6)} {mp.nstr(W1,5)} | {[round(float(ob[b]*Qb[b]-d),1) for b in range(1,bmax+1)]}")

# monotonicity check of gstar in Q on a grid for k=12..17 (Q from Q_bmax to 1e4)
P("\nmonotonicity of gstar(Q) on grid Q in [Q_bmax, 1e4] (should be decreasing):")
for k in range(12, 18):
    n, d, r, w, j, mp_ = params(k)
    bmax = (mp_ - 1) // (2 * j)
    o1 = mp_ - 2 * j
    N2 = Fr(n * (n - 1), 2); Pm = (2**bmax) * Fr((n // 2**bmax) * (n // 2**bmax - 1), 2); Qm = float(N2 / Pm)
    if o1 * Qm <= d:
        P(f"k={k}: o1*Qm <= d, case B gives no failure"); continue
    def gs(Q):
        Q = mp.mpf(Q)
        return 2 * (mp_ * mp.log(4) + o1 * mp.log(Q * n / (n - 1)) - mp.log(mp.mpf(d) / (mp_ + 1))) / (o1 * Q - d)
    grid = [Qm * (1.01**i) for i in range(0, 700)]
    vals = [gs(Q) for Q in grid]
    dec = all(vals[i + 1] < vals[i] for i in range(len(vals) - 1))
    P(f"k={k}: decreasing={dec}, gstar(Qm)={mp.nstr(vals[0],6)}, gstar(1e4-ish)={mp.nstr(vals[-1],6)}")

open(__file__.replace('.py', '_out.txt'), 'w').write("\n".join(out) + "\n")
