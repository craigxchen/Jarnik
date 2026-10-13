"""Judge checks J2/J3: where the leading-order (tangent-cone) criterion can work.

For the slab box of round10 Thm 2.1 with box size T (T = n as in round 10, or T = 1):
  m(n,T) = largest m with 8 binom(m+n-2, n-1) <= (2T+1)^(n-1)          (EXACT big-integer binomials)
  w      = T sqrt(8 (n-1) ln 4r)
  D      = 2 (n-1) T
J2 (Chebyshev): for c in V(S,m), |Phi_c(u)| <= |c|_1 * 2 (R/2)^m with R = max_S |<k',u>| <= D/2
   (|u|_inf = 1), so kappa(u) <= 2^(1-m). Criterion (2.1) of theory.md needs
   kappa > (D/2)|u|_inf/(m+1), and |u|_inf >= e^(-W/2) for distinct lattice points, so (2.1) can
   hold only if  W > W_TC := 2(m-1) ln 2 + 2 ln((D/2)/(m+1)).
   Printed: W_TC against W_min(n) = 3 n ln n - 14 n (two_thirds lower bound, C = 1).
J3 (T = 1 re-parametrisation of round10 Thm 2.1/4.1):
   * m >= 4w threshold (r = n, n ln n, n^2) and n - 1 >= 300 ln 4r;
   * Theorem 4.1 comparison: (m/2 - w) W <= 2 ln|c|_1 + 2m ln(D C/2) - 2 ln m!, with
     |c|_1 <= U H = U^3 (H = (U T^(m-1))^2 = U^2 at T = 1), U <= 3^(n-1).  W_1 = that bound / (m/2 - w).
   * room for (2.1) at the actual scale |u|_inf <= C e^(-W/4):
     log room = ln(2^(1-m)) - ln((D/2) e^(-W/4)/(m+1)) at W = W_min(n), C = 1.
"""
import math, sys, functools


@functools.lru_cache(maxsize=None)
def m_exact(n, T):
    rhs = (2 * T + 1) ** (n - 1)
    # m upper bound: binom grows; search
    lo, hi = 0, 4 * n * T + 10
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if 8 * math.comb(mid + n - 2, n - 1) <= rhs:
            lo = mid
        else:
            hi = mid - 1
    return lo


def lnfact(m):
    return math.lgamma(m + 1)


def rowT1(n, r):
    T = 1
    m = m_exact(n, T)
    w = T * math.sqrt(8 * (n - 1) * math.log(4 * r))
    D = 2 * (n - 1) * T
    lnU = (n - 1) * math.log(3)
    lnc1 = 3 * lnU
    num = 2 * lnc1 + 2 * m * math.log(D / 2) - 2 * lnfact(m)
    W1 = num / (m / 2 - w) if m / 2 > w else float('inf')
    Wmin = 3 * n * math.log(n) - 14 * n
    WTC = 2 * (m - 1) * math.log(2) + 2 * math.log((D / 2) / (m + 1))
    room = (1 - m) * math.log(2) - (math.log(D / 2) - Wmin / 4 - math.log(m + 1))
    return m, w, W1, Wmin, WTC, room


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "Tn"):
        print("J2: T = n (round 10 parameters), r = n")
        print("n, m(exact), m/(n(n-1)), W_TC = 2(m-1)ln2 + 2ln((D/2)/(m+1)), W_min = 3n ln n - 14n, W_TC/W_min")
        for n in [64, 256, 1024, 2137, 3000]:
            T = n
            m = m_exact(n, T)
            D = 2 * (n - 1) * T
            WTC = 2 * (m - 1) * math.log(2) + 2 * math.log((D / 2) / (m + 1))
            Wmin = 3 * n * math.log(n) - 14 * n
            print(f"{n}, {m}, {m/(n*(n-1)):.4f}, {WTC:.6g}, {Wmin:.6g}, {WTC/Wmin:.1f}")
        print()
    if which in ("all", "T1"):
        print("J3: T = 1, exact m; first n (scanning upward) from which m >= 4w holds for every n up to the scan end")
        for label, rf in [("r=n", lambda n: n), ("r=n ln n", lambda n: int(n * math.log(n)) + 1), ("r=n^2", lambda n: n * n)]:
            first_ok = None
            last_fail = None
            lo, hi = 1000, 9000
            for n in range(lo, hi + 1, 1):
                r = rf(n)
                m = m_exact(n, 1)
                w = math.sqrt(8 * (n - 1) * math.log(4 * r))
                if m >= 4 * w:
                    if first_ok is None:
                        first_ok = n
                else:
                    last_fail = n
            print(f"  {label}: m>=4w first at n={first_ok}; last failure in [{lo},{hi}] at n={last_fail}; "
                  f"300 ln 4r threshold: n-1 >= {300*math.log(4*rf(last_fail+1 if last_fail else lo)):.0f} at n={last_fail+1 if last_fail else lo}")
        print()
        print("J3 table (r = n): n, m(T=1), m/n, w, m/(4w), W_1 (Thm 4.1 at T=1, C=1), W_min, W_TC, log room at W_min")
        for n in [1000, 2000, 3000, 3500, 4000, 5000, 6000, 8000]:
            m, w, W1, Wmin, WTC, room = rowT1(n, n)
            print(f"{n}, {m}, {m/n:.4f}, {w:.1f}, {m/(4*w):.3f}, {W1:.2f}, {Wmin:.6g}, {WTC:.6g}, {room:.6g}")
