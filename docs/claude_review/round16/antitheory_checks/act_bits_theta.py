"""AC_T at threshold theta: Q delta-small on a set of measure >= theta.  theta = 1/2 (steelman),
theta0 = 1 - beta0, beta0 = 2r(4r)^{-3/2} (profile-free lower bound for mu(A): sub-Gaussian, sum c_x^2 <= 4d),
thetaS = 1 - betaS, betaS = 2r exp(-3 w^2/(8n)) (Sylvester sign columns: sum c_x^2 = 2n).
Bit strata with |X_i| = n/2: exceptional iff rho = 2cq <= m'-1 with q*tail(c) <= 1 - theta.
Gamma test at W in {W_min, 3n ln n, 4 n ln n} for the y + eta v point (see act_strata.py)."""
import sys
from common import *
from act_strata import best_c, gamma_point

def run(kmin, kmax, gamma_kmax):
    for k in range(kmin, kmax + 1):
        n = 2 ** k; r = n - 1
        d, w, j, m = wjm(n, r)
        beta0 = 2 * r * (4 * mpf(r)) ** (-1.5)
        betaS = 2 * r * exp(-3 * w ** 2 / (8 * n))
        for name, beta in [('half', mpf(1) / 2), ('theta0', beta0), ('thetaS', betaS)]:
            best = None; fails = []
            for q in range(1, k):
                c = best_c(n // 2, beta / q)
                rho = 2 * c * q
                if rho > m - 1:
                    continue
                o = m - rho; K = 2 ** q
                nec = K * o > d
                tag = ''
                if nec and k <= gamma_kmax:
                    for W, wn in [(Wmin(n), 'Wmin'), (3 * n * log(n), '3nlnn'), (4 * n * log(n), '4nlnn')]:
                        le = ln_inv_eps0(n, m, W)
                        eta = exp(-(le + m * log(4)) / o)
                        if eta > mpf(1) / (2 * (K - 1)):
                            continue
                        if gamma_point(n, q, eta) <= W / (2 * d):
                            tag = 'FAIL@' + wn; break
                fails.append((q, c, rho, o, nec, tag))
            qmax = max([f[0] for f in fails], default=0)
            bad = [(f[0], f[5]) for f in fails if f[4]]
            print(k, name, 'j=%d c/j(q=4)=%s' % (j, mp.nstr(mpf(best_c(n // 2, beta / 4)) / j, 3)), 'q_max=', qmax, 'K*o>d at q in', [b[0] for b in bad], 'Gamma-fail:', [b for b in bad if b[1]])
        sys.stdout.flush()

if __name__ == '__main__':
    run(12, 22, 22)
    run(23, 60, 0)
