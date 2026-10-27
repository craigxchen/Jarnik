"""What CIRCLE(q) 'yes' would require of Psi_A(C_1) at regime sdc_q, through kappa_V (EVIDENCE-level bookkeeping; the
formulas are exact, the moment model in (c) is heuristic).  C = 1, c_0 = 0.184 (H_light evidence value), tau = tau_D.
 (a) threshold T = eps_0/(2 sqrt(m'+1) e^{pi^2/(4q)}),  eps_0 = C d e^{-W/4}/(m'+1),  W = 2 q tau_D,  d = q;
 (b) kappa_V(hat u) <= T  forces  ln Psi_A(C_1) <= 2 (ln T + m' ln d)                     [by (2.1)]
     and is implied by  ln Psi_A(C_1) <= 2 (ln T + m' ln d) - ln U,  U <= 3^d             [by (1.1)];
 (c) heuristic size of E_A l^{2m'} on C_1: Gaussian model with variance q/3 (cube value; regime A has density >= 1/2):
     ln E ~ m' ln(q/3) + ln (2m'-1)!!;  required relative residual  -ln rel_need / m'.
Also ln of the cube's relative residual on C_1 in the Gaussian model:  ln(m!/(2m-1)!!) ~ -m ln 2.
"""
import math
import mpmath as mp


def tau_D(q, c0=0.184, C=1.0):
    r = 2 * q; n = q + 1
    return max((2 / c0) * (2 * math.log(4 * r + 5) + math.log(1000 / C)), 4 * math.log(4 * r + 5) + 8,
               2 * math.log(1000 * q * n * n / C))


for q in (3251, 5479, 10007, 100003, 10 ** 6 + 3):
    d = q; r = 2 * q
    w = math.sqrt(8 * d * math.log(4 * r)); j = math.floor(w / 2); mp_ = math.ceil(4 * w)
    tD = tau_D(q); W = 2 * q * tD
    ln_eps0 = math.log(d / (mp_ + 1)) - W / 4
    lnT = ln_eps0 - math.log(2 * math.sqrt(mp_ + 1)) - math.pi ** 2 / (4 * q)
    nec = 2 * (lnT + mp_ * math.log(d)); suf = nec - d * math.log(3)
    lnE = mp_ * math.log(q / 3) + float(mp.log(mp.fac2(2 * mp_ - 1)))
    need_rate = (lnE - nec) / mp_
    cube_rate = -float(mp.log(mp.factorial(mp_) / mp.fac2(2 * mp_ - 1))) / mp_
    print(f"q={q}: w={w:.2f} j={j} m'={mp_} m'/q={mp_/q:.3f} tau_D={tD:.1f} W=2q tau_D={W:.4e} ln eps0={ln_eps0:.4e} ln T={lnT:.4e}")
    print(f"     kappa_V-cap needs ln Psi_A(C_1) <= {nec:.4e} (necessary), <= {suf:.4e} (sufficient);"
          f" heuristic ln E_A l^2m' ~ {lnE:.4e}; required decay -ln(rel)/m' >= {need_rate:.1f} per degree;"
          f" cube (Gaussian model) decay = {cube_rate:.3f} per degree")
