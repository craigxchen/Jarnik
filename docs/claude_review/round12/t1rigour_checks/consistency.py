"""Scale-consistency of the T = 1, m' = ceil(4w) reduction (t1rigour.md §5.3, Prop 5.4).
For a cheap family with s = sum(ell_i) and exceptional set of size e, the separation lemma gives radius
  sigma = (1/2) exp(-Q abar/2),  Q = 2^s n(n-1)/((n-e)(n-e-2^s)),  abar = W/(2(n-1)) + 2 ln C',
and the cost satisfies d_F >= 2 j s + 2 e. The transversal upper bound kappa <= 3^m' eta^(m'-d_F) is compatible
with the requirement kappa > eps0 = C d e^(-W/4)/(m'+1) at eta = sigma iff
  G := m' ln 3 - (m' - 2js - 2e)(ln 2 + Q abar/2) - ln(C d/(m'+1)) + W/4 > 0.
We scan all admissible (s, e) and report the first n (r = n, C = 1) from which min G > 0, for W = W_min = 3n ln n - 14 n
and for the W -> infinity limit (coefficient of W: 1/4 - (m'-2js-2e) Q/(4(n-1)) > 0).
"""
import math

def check(n, W=None, C=1.0):
    d = n - 1
    w = math.sqrt(8 * d * math.log(4 * n))
    mp = math.ceil(4 * w)
    j = math.floor(w / 2)
    b = (mp - 1) // (2 * j)
    Cp = max(C, 1.0)
    worst = float('inf')
    for s in range(0, b + 1):
        emax = (mp - 1) // 2 - j * s
        for e in range(0, max(emax, 0) + 1):
            if n - e - 2 ** s <= 0:
                return -1.0
            Q = 2 ** s * n * (n - 1) / ((n - e) * (n - e - 2 ** s))
            order = max(mp - 2 * j * s - 2 * e, 0)
            if W is None:
                G = 0.25 - order * Q / (4 * (n - 1))
            else:
                abar = W / (2 * (n - 1)) + 2 * math.log(Cp)
                G = mp * math.log(3) - order * (math.log(2) + Q * abar / 2) - math.log(C * d / (mp + 1)) + W / 4
            worst = min(worst, G)
    return worst

def first(fn, lo=100, hi=10**7):
    # last failure on a scan
    last = None
    n = lo
    while n <= hi:
        if fn(n) <= 0:
            last = n
        n = n + 1 if n < 20000 else int(n * 1.001) + 1
    return last

print("last n (r=n, C=1) with min G <= 0, W -> infinity:", first(lambda n: check(n, None)))
print("last n (r=n, C=1) with min G <= 0, W = 3n ln n - 14n:", first(lambda n: check(n, 3 * n * math.log(n) - 14 * n)))
for n in (3000, 5000, 6000, 10000, 10**5):
    print(n, "G_inf =", round(check(n, None), 4), " G(W_min) =", round(check(n, 3 * n * math.log(n) - 14 * n), 2))
