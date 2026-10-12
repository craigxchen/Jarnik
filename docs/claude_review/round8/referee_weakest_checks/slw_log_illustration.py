"""Referee check R2 (illustration): systems that SLW^log(theta, g, kappa) must cover, violating it.

Vertical templates with D = 2s columns b = 1..D, s odd-power conditions, k = s integer kernel forms
(LLL basis).  gamma -> (2s+1)/(4s) > 1/2.  For each system we print k, n (= number of primes), the
SLW^log threshold theta*log(n+2) at the largest admissible theta = 1/log 5.5, the plain exponent gamma,
and the "effective" exponent including c = 1/2 and kappa = 1 (all coefficient factors, B_0 included).
The asymptotic refutation (t -> infinity with s growing slowly) is Prop. R2 of the report; these
rows only show that the covered systems already violate the inequality at moderate t.
"""
import mpmath
from slw_template_counterexample import analyse, kernel_basis, modulus_Q, int_rank

mpmath.mp.dps = 500


if __name__ == '__main__':
    theta = 1 / mpmath.log(5.5)
    for s in (3, 4, 5, 6):
        D = 2 * s
        b = list(range(1, D + 1))
        sig = kernel_basis(b, s)
        assert len(sig) == s and int_rank(sig) == s
        Q = modulus_Q(b)
        maxc = max(abs(x) for v in sig for x in v)
        for u in (101, 10007):
            t = Q * u
            r = analyse(b, sig, t)
            assert r['rank'] == s
            thr = theta * mpmath.log(r['n'] + 2)
            geff = (mpmath.log(1 / r['lam']) - (r['n'] + 1) * mpmath.log(2) - r['coef']) / r['logH']
            print(f"s={s} D={D} k={s} max|sigma|={maxc:<4d} log t={float(r['logt']):6.2f} n={r['n']:3d} "
                  f"theta*log(n+2)={float(thr):5.2f} covered(k>=thr)={s >= thr} "
                  f"gamma={float(r['gamma']):.4f} (limit {(2*s+1)/(4*s):.4f}) gamma_eff(c=1/2,kappa=1)={float(geff):.4f}")
