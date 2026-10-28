#!/usr/bin/env python3
"""Referee: uncut cube {-1,0,1}^q on the circle u_x = cos(2 pi x/q), at large (q, m).
Psi_cube = (m!)^2 [t^m] prod_x (1 + (2/3) u_x^2 t + u_x^4 t^2/18)   (reduction Thm 5.1(a); verified exactly above for q=7)
E_cube l^{2m} = (2m)! [s^{2m}] prod_x (1 + 2 cosh(u_x s))/3.
Both via log/exp power series (exact power sums of cos^{2n}, valid for 2n < 2q, i.e. n < q), mpmath 40 digits.
Prints ln Psi_cube, ln E, -ln(rel)/m, and the Gaussian-model ln E used in calib §7.1."""
import sys, time
import mpmath as mp
import os
mp.mp.dps = int(os.environ.get("DPS","40"))

def powsum(q, n):
    # sum_x cos(2 pi x/q)^{2n} = q binom(2n,n)/4^n  for n < q (q odd prime)
    assert n < q
    return q * mp.binomial(2*n, n) / mp.mpf(4)**n

def exp_series(s, N):
    # e = exp(S), S = sum_{k>=1} s[k] z^k ; n e_n = sum_{k=1}^n k s_k e_{n-k}
    e = [mp.mpf(0)]*(N+1); e[0] = mp.mpf(1)
    for n in range(1, N+1):
        acc = mp.mpf(0)
        for k in range(1, n+1):
            if s[k] != 0:
                acc += k * s[k] * e[n-k]
        e[n] = acc / n
    return e

def run(q, m):
    t0 = time.time()
    # Psi: log(1 + a t + b t^2) = log(1-al t) + log(1-be t), al,be = roots-reciprocals: al+be = -a, al be = b
    # per unit u^2: al0, be0 with a=2/3, b=1/18
    disc = mp.sqrt(mp.mpf(2)/9)
    al0 = (-mp.mpf(2)/3 + disc)/2; be0 = (-mp.mpf(2)/3 - disc)/2
    # check: al0+be0 = -2/3, al0*be0 = (4/9 - 2/9)/4 = 1/18
    s = [mp.mpf(0)]*(m+1)
    for n in range(1, m+1):
        s[n] = -(al0**n + be0**n)/n * powsum(q, n)
    e = exp_series(s, m)
    Psi = mp.factorial(m)**2 * e[m]
    # E l^{2m}: log((1+2cosh y)/3) = sum_k c_k y^{2k}; get c_k by series log of 1 + sum_{k>=1} 2 y^{2k}/(3 (2k)!)
    g = [mp.mpf(0)]*(m+1)
    for k in range(1, m+1):
        g[k] = 2/(3*mp.factorial(2*k))
    # log(1+G): l_n via n l_n = n g_n - sum_{k=1}^{n-1} k l_k g_{n-k}
    lg = [mp.mpf(0)]*(m+1)
    for n in range(1, m+1):
        acc = n*g[n]
        for k in range(1, n):
            acc -= k*lg[k]*g[n-k]
        lg[n] = acc/n
    s2 = [mp.mpf(0)]*(m+1)
    for k in range(1, m+1):
        s2[k] = lg[k]*powsum(q, k)
    e2 = exp_series(s2, m)
    E = mp.factorial(2*m) * e2[m]
    rel = Psi/E
    gauss = mp.log(mp.factorial(2*m)/(mp.mpf(2)**m*mp.factorial(m))) + m*mp.log(mp.mpf(q)/3)
    print(f'q={q} m={m}: ln Psi_cube={mp.nstr(mp.log(Psi),10)}  ln E_cube l^2m={mp.nstr(mp.log(E),10)} '
          f'(Gaussian model {mp.nstr(gauss,10)})  -ln(rel)/m={mp.nstr(-mp.log(rel)/m,8)}  m/q={m/q:.3f}  ({time.time()-t0:.0f}s)', flush=True)

if __name__ == '__main__':
    for a in sys.argv[1:]:
        q, m = map(int, a.split(':'))
        run(q, m)
