# Own recomputation: K.gamma-hat and (K+D_PQ).gamma-hat on Mbar_{0,k+2} (marks 1..k,P,Q)
# gamma-hat: D_PQ -> 1/4, D_{P u A} (nonempty proper A of [k]) -> 2^-k, all other boundary divisors -> 0
# K = sum_S (s(n-s)/(n-1) - 2) D_S  (Keel-McKernan / Pandharipande formula), n=k+2
from fractions import Fraction as F
from math import comb
for k in range(3,10):
    n=k+2
    c=lambda s: F(s*(n-s),n-1)-2
    K = F(1,4)*c(2) + sum(comb(k,a)*c(a+1) for a in range(1,k))*F(1,2**k)
    KD = K + F(1,4)
    cnt = 2**(k-3)*(5-k)-1 if k>=3 else None
    # direct count exponent per w: W'(1-(k-1)(1/4+delta/2)), delta=1/(2^k-2), W'=(2^(k-1)-1) w
    d=F(1,2**k-2); direct=(2**(k-1)-1)*(1-(k-1)*(F(1,4)+d/2))
    # Vojta margin kappa at delta_k and D=0 margin
    kap=(k-1)*(1+2*d)/4-1; kap0=(k-2)*(1+2*d)/4-1
    print(k,'K.g=',K,' (K+D).g=',KD,' count=',cnt,' direct=',direct,' -(K+D).g*2^(k-1)=',-KD*2**(k-1),' kappa(delta_k)=',kap,' D=0 margin=',kap0)
