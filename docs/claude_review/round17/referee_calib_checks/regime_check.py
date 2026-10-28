#!/usr/bin/env python3
"""Referee: re-derive calib §7.1 bookkeeping (C=1, c0=0.184, tau=tau_D of barrier Lemma 3.4, d=q, r=2q, n=q+1)."""
from math import log, sqrt, ceil, floor, pi, lgamma
from sympy import isprime
def row(q, c0=0.184, C=1.0, tau=None):
    n = q+1; r = 2*q; d = q
    w = sqrt(8*d*log(4*r)); j = floor(w/2); m = ceil(4*w)
    tauD = max((2/c0)*(2*log(4*r+5)+log(1000/C)), 4*log(4*r+5)+8, 2*log(1000*q*n*n/C))
    t = tauD if tau is None else tau
    W = 2*q*t
    lneps0 = log(C*d) - W/4 - log(m+1)
    lnT = lneps0 - log(2) - 0.5*log(m+1) - pi**2/(4*q)
    need = 2*(lnT + m*log(d))
    gauss = lgamma(2*m+1) - m*log(2) - lgamma(m+1) + m*log(q/3)
    return dict(q=q, j=j, m=m, mq=m/q, tauD=tauD, W=W, lnT=lnT, need=need, gauss=gauss, rate=(gauss-need)/m,
                regime=(d >= 316*log(4*r)), m2j1=m/(2*j+1))
for q in [3251, 5479, 10007, 100003, 1000003]:
    R = row(q); print({k: (round(v,4) if isinstance(v,float) else v) for k,v in R.items()})
R = row(3251, tau=2*log(4*2*3251)); print('tau floor 2ln4r:', round(2*log(4*2*3251),2), 'rate', round(R['rate'],2))
# cube-model prediction at q=3251 using exact cube rel (ln Psi_cube=27188.36, ln E_cube=28949.37)
R = row(3251); lnk = 0.5*27188.36 - R['m']*log(3251)
print('cube-model ln k_lo at q=3251:', round(lnk,1), ' ln T:', round(R['lnT'],1), ' gap:', round(lnk-R['lnT'],1))
# first regime prime q=3 mod 4
q=3
while True:
    q+=1
    if q%4==3 and isprime(q) and q >= 316*log(8*q): print('first regime q =', q); break
