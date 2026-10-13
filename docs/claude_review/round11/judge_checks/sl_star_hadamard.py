"""Judge check J1: Hypothesis SL* fails for Sylvester-Hadamard profiles.

Part 1 (exact): for every k in 3..7, every pair of distinct nonzero a, b in F_2^k and every
class c in F_2^2 of (a.x, b.x), there are integers (alpha, beta, gamma, delta) with
    1[a.x=0] = alpha + beta 1[b.x=0] + gamma 1[(a+b).x=0] + delta 1[(a.x,b.x)=c]   for all x,
with |delta| = 2. (So a column is a Z-combination of 1, two other columns, and the indicator of
n/4 points; choosing c != class of the base point keeps the base point out of the support.)

Part 2 (exact integers / proved bound): cheapness of that family in the regime slab box
(T = n, w = n sqrt(8 (n-1) ln 4r), r = n - 1):
    cost = 2 * (2 floor(w/2)) + (n/4) * 2T   <   m,
using the PROVED lower bound m >= 0.7 n (n-1) (round10 Thm 2.1(b)), and also with m computed
from log-gamma at 60 digits (labelled). Also prints the regime conditions.
"""
import itertools, math
from fractions import Fraction
import mpmath

mpmath.mp.dps = 60


def dot(a, x):
    return sum(ai * xi for ai, xi in zip(a, x)) % 2


def part1():
    out = []
    for k in range(3, 8):
        pts = list(itertools.product([0, 1], repeat=k))
        nz = [v for v in pts if any(v)]
        bad = 0
        checked = 0
        for a in nz:
            for b in nz:
                if a == b:
                    continue
                ab = tuple((ai + bi) % 2 for ai, bi in zip(a, b))
                for cls in itertools.product([0, 1], repeat=2):
                    # solve on the four classes: s = al + be t + ga u + de [class]
                    # classes (a.x, b.x) -> (s, t, u)
                    rows = []
                    for (ax, bx) in itertools.product([0, 1], repeat=2):
                        s = 1 if ax == 0 else 0
                        t = 1 if bx == 0 else 0
                        u = 1 if (ax + bx) % 2 == 0 else 0
                        ind = 1 if (ax, bx) == cls else 0
                        rows.append(([1, t, u, ind], s))
                    # 4x4 linear system over Q
                    M = [[Fraction(v) for v in r[0]] + [Fraction(r[1])] for r in rows]
                    # gaussian elimination
                    n = 4
                    for col in range(n):
                        piv = next(i for i in range(col, n) if M[i][col] != 0)
                        M[col], M[piv] = M[piv], M[col]
                        for i in range(n):
                            if i != col and M[i][col] != 0:
                                f = M[i][col] / M[col][col]
                                M[i] = [x - f * y for x, y in zip(M[i], M[col])]
                    sol = [M[i][n] / M[i][i] for i in range(n)]
                    # verify pointwise on F_2^k
                    ok = all(sol[i].denominator == 1 for i in range(4)) and abs(sol[3]) == 2
                    for x in pts:
                        s = 1 if dot(a, x) == 0 else 0
                        t = 1 if dot(b, x) == 0 else 0
                        u = 1 if dot(ab, x) == 0 else 0
                        ind = 1 if (dot(a, x), dot(b, x)) == cls else 0
                        if sol[0] + sol[1] * t + sol[2] * u + sol[3] * ind != s:
                            ok = False
                    checked += 1
                    if not ok:
                        bad += 1
                    if k == 3 and a == nz[0] and b == nz[1]:
                        out.append(f"  example k=3 a={a} b={b} class={cls}: (alpha,beta,gamma,delta) = {tuple(str(v) for v in sol)}")
        out.append(f"k={k}: n={2**k}, checked {checked} (a,b,class) triples, failures {bad}")
    return out


def log_binom(N, K):
    return mpmath.loggamma(N + 1) - mpmath.loggamma(K + 1) - mpmath.loggamma(N - K + 1)


def m_exact_lgamma(n, T):
    """largest m with 8 binom(m+n-2, n-1) <= (2T+1)^(n-1), via 60-digit log-gamma."""
    rhs = (n - 1) * mpmath.log(2 * T + 1) - mpmath.log(8)
    lo, hi = 1, 4 * n * T + 10
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if log_binom(mid + n - 2, n - 1) <= rhs:
            lo = mid
        else:
            hi = mid - 1
    # margin: distance of both neighbours from the threshold (to show the decision is robust)
    marg_in = rhs - log_binom(lo + n - 2, n - 1)
    marg_out = log_binom(lo + 1 + n - 2, n - 1) - rhs
    return lo, float(marg_in), float(marg_out)


def part2():
    out = []
    out.append("k, n, w, 2*2floor(w/2), (n/4)*2T, cost, 0.7n(n-1) [proved m lower bound], cost<0.7n(n-1)?, "
               "m(lgamma), cost<m?, n-1>=300ln4r?, m>=4w (with proved bound)?")
    for k in range(11, 21):
        n = 2 ** k
        r = n - 1
        T = n
        w = T * math.sqrt(8 * (n - 1) * math.log(4 * r))
        col = 2 * (2 * math.floor(w / 2))
        coords = (n // 4) * 2 * T
        cost = col + coords
        mlow = 0.7 * n * (n - 1)
        regime = (n - 1) >= 300 * math.log(4 * r)
        m4w = mlow >= 4 * w
        if k <= 16:
            m, mi, mo = m_exact_lgamma(n, T)
            mstr = f"{m} (margins {mi:.3g},{mo:.3g})"
            cm = cost < m
        else:
            mstr, cm = "-", "-"
        out.append(f"{k}, {n}, {w:.6g}, {col}, {coords}, {cost:.6g}, {mlow:.6g}, {cost < mlow}, {mstr}, {cm}, {regime}, {m4w}")
    return out


if __name__ == "__main__":
    print("Part 1: column identity on F_2^k (exact)")
    for line in part1():
        print(line)
    print()
    print("Part 2: cheapness in the regime slab box (T = n, r = n - 1)")
    for line in part2():
        print(line)
