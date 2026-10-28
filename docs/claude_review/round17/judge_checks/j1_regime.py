# J1: regime arithmetic for CIRCLE(q) with the now-PROVED H_light constants (C = 1).
# 60-digit mpmath; all quantities are explicit closed forms, so this is exact up to printed rounding.
from mpmath import mp, mpf, log, sqrt, ceil, floor, pi, exp
mp.dps = 60
C = mpf(1)
def regime(q):
    q = mpf(q); d = q; n = q + 1; r = 2*q
    reg = d >= 316*log(4*r)
    w = sqrt(8*d*log(4*r)); j = int(floor(w/2)); mp_ = int(ceil(4*w))
    return dict(q=q, d=d, n=n, r=r, reg=reg, w=w, j=j, m=mp_)
def tauD(R, c0):
    q, n, r = R['q'], R['n'], R['r']
    t1 = (2/c0)*(2*log(4*r+5) + log(1000/C)); t2 = 4*log(4*r+5) + 8; t3 = 2*log(1000*q*n**2/C)
    return max(t1, t2, t3), (t1, t2, t3)
def at(R, tau):
    q, d, m = R['q'], R['d'], R['m']
    W = 2*q*tau
    lneps0 = log(C*d) - W/4 - log(m+1)
    lnT = lneps0 - log(2) - log(m+1)/2 - pi**2/(4*q)
    budget = q*tau/2 - m*log(2*q) - log(C*q/(m+1))      # survive Cor 3.2: ln lambda_q must be below this
    needPsi = 2*(lnT + m*log(d))                         # calib (2.1): kappa_V <= T forces ln Psi_A(C_1) <= this
    lnHp = m*W/40
    return dict(W=W, lnT=lnT, budget=budget, budget_q=budget/q, needPsi=needPsi, lnHp=lnHp)
for qq in [3251, 5479]:
    R = regime(qq)
    q = R['q']
    LB = min(max(q/3, q/2 - 11 - mpf(3)/2*sqrt(q)), (q-3)/2)
    cands = [('1/6 (elementary, all q)', mpf(1)/6), ('LB(q)/r (Weil, Thm 2.8)', LB/R['r'])]
    if qq == 3251: cands.append(('exact min 1571/r', mpf(1571)/R['r']))
    if qq == 5479: cands.append(('exact min 2667/r', mpf(2667)/R['r']))
    cands.append(('0.184 (R16 evidence value, superseded)', mpf('0.184')))
    print(f"q={qq}: regime={R['reg']} w={mp.nstr(R['w'],8)} j={R['j']} m'={R['m']} m'<q={R['m']<qq} LB(q)={mp.nstr(LB,8)}")
    for name, c0 in cands:
        tD, parts = tauD(R, c0); A = at(R, tD)
        print(f"  c0={mp.nstr(c0,6):>8} [{name}]: tau_D={mp.nstr(tD,7)} (terms {mp.nstr(parts[0],6)}, {mp.nstr(parts[1],6)}, {mp.nstr(parts[2],6)})")
        print(f"      ln T* = {mp.nstr(A['lnT'],8)}; survival budget ln|c|_1 < {mp.nstr(A['budget'],8)} = {mp.nstr(A['budget_q'],6)} q;"
              f" kappa_V-cap needs ln Psi_A(C_1) <= {mp.nstr(A['needPsi'],8)}; ln H' = {mp.nstr(A['lnHp'],6)}")
    # height window of survive Thm 3.4 at the proved c0: [ln H_*, ln H'] with ln H_* = budget - ln|A|, |A| <= 3^q
    tD,_ = tauD(R, LB/R['r']); A = at(R, tD)
    print(f"  height window at proved c0 (Weil): ln H in [{mp.nstr(A['budget']-q*log(3),7)} (using |A|<=3^q), {mp.nstr(A['lnHp'],7)}]")
    # c0-free floor (third term), as in cap.md
    t3 = tauD(R, mpf(1))[1][2]; A = at(R, t3)
    print(f"  c0-free floor tau={mp.nstr(t3,7)}: ln T*={mp.nstr(A['lnT'],8)}, budget={mp.nstr(A['budget_q'],6)} q")
# needed average relative-residual decay per degree for a kappa_V-cap (calib §7.1 quantity), q = 3251,
# using calib's Gaussian-model ln E l^{2m'} = 2.94e4 (calib regime_need.txt) and the referee's exact cube value 28949.
R = regime(3251); q = R['q']
for name, c0 in [('Weil LB', (q/2 - 11 - mpf(3)/2*sqrt(q))/R['r']), ('exact', mpf(1571)/R['r']), ('0.184', mpf('0.184'))]:
    A = at(R, tauD(R, c0)[0])
    for lnE in [mpf(29437), mpf(28949)]:
        print(f"decay needed (c0 {name}, lnE {int(lnE)}): {mp.nstr((lnE - A['needPsi'])/R['m'], 6)} per degree")
