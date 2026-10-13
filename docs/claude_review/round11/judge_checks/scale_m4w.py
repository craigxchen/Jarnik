"""Judge check J3b: T = 1 with the smallest LP-admissible order m' = ceil(4w) (instead of the maximal m).

Columns (r = n, C = 1):
  w = sqrt(8 (n-1) ln 4n);  m' = ceil(4w);  m_max(n, T=1) from 60-digit log-gamma (largest m with
  8 binom(m+n-2, n-1) <= 3^(n-1)); need m' <= m_max for V(A,m') to contain V(A,m_max) (it does: V(A,m) grows as m
  decreases) and for the Siegel count (E <= U/4 holds a fortiori).
  W_min = 3 n ln n - 14 n (two_thirds lower bound, C = 1).
  W_1   = Thm 4.1 threshold at T = 1 with order m': [6 (n-1) ln 3 + 2 m' ln(n-1) - 2 ln m'!] / (m'/2 - w).
  cols  = number of e_p = 1 columns that are cheap together: largest k with k * w < m'  (h_A(a_q) <= 2 floor(w/2) <= w).
  rho0  = e^(-W_min/(4 m')): heuristic radius of TC_eps around a cheap span of small cost (polynomial factors dropped).
  eta_t = e^(-W_min/(4(n-1))) / sqrt(n): pair-inequality separation for a Plotkin-tight pair (excess W/(2(n-1))).
  ok    = eta_t > rho0 (heuristic comparison; constants not optimised).
Everything except m_max is elementary floating point; the comparison is a heuristic scale count (README §5.3).
"""
import math
import mpmath

mpmath.mp.dps = 60


def m_max_T1(n):
    rhs = (n - 1) * mpmath.log(3) - mpmath.log(8)
    lo, hi = 1, 2 * n
    while lo < hi:
        mid = (lo + hi + 1) // 2
        lb = mpmath.loggamma(mid + n - 1) - mpmath.loggamma(n) - mpmath.loggamma(mid)
        if lb <= rhs:
            lo = mid
        else:
            hi = mid - 1
    return lo


print("n, w, m'=ceil(4w), m_max(T=1), W_min, W_1(m'), cheap cols, ln rho0, ln eta_t, eta_t>rho0")
for n in [3000, 5000, 10000, 20000, 50000, 100000, 1000000]:
    w = math.sqrt(8 * (n - 1) * math.log(4 * n))
    mp = math.ceil(4 * w)
    mm = m_max_T1(n)
    Wmin = 3 * n * math.log(n) - 14 * n
    W1 = (6 * (n - 1) * math.log(3) + 2 * mp * math.log(n - 1) - 2 * math.lgamma(mp + 1)) / (mp / 2 - w)
    cols = max(k for k in range(0, 100) if k * w < mp)
    lnrho = -Wmin / (4 * mp)
    lneta = -Wmin / (4 * (n - 1)) - 0.5 * math.log(n)
    print(f"{n}, {w:.1f}, {mp}, {mm}, {Wmin:.6g}, {W1:.1f}, {cols}, {lnrho:.2f}, {lneta:.2f}, {lneta > lnrho}")
