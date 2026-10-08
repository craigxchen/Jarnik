"""Exact minimum of the level deficit D = sum_tau (S_tau - M/2)^2 over
allocations of M points to consecutive levels with every level occupancy
<= c (first/last level nonempty), compared with the bound
D >= (M-4c)^3/(12c)  (used in Theorem A of cotangent.md)."""
from fractions import Fraction
from functools import lru_cache

def min_deficit(M, c):
    # DP over partial sums S (points already placed below the threshold)
    INF = float('inf')
    @lru_cache(maxsize=None)
    def f(S):
        # S points placed in levels below current threshold (S>=1); choose next
        # level size n (1..c); if S+n == M we stop (that level is the top level).
        best = INF
        for n in range(1, c + 1):
            T = S + n
            if T == M:
                best = min(best, 0)
            elif T < M:
                best = min(best, (Fraction(T) - Fraction(M, 2)) ** 2 + f(T))
        return best
    best = INF
    for n0 in range(1, min(c, M) + 1):
        if n0 == M:
            best = 0
        else:
            best = min(best, (Fraction(n0) - Fraction(M, 2)) ** 2 + f(n0))
    return best

bad = 0
worst = None
for M in range(2, 90):
    for c in range(1, M + 1):
        D = min_deficit(M, c)
        bound = Fraction(max(M - 4 * c, 0) ** 3, 12 * c)
        if D < bound:
            bad += 1
        if M >= 8 * c:
            r = D / Fraction(M ** 3, 12 * c)
            worst = r if worst is None or r < worst else worst
print('violations of D >= (M-4c)^3/(12c):', bad)
print('min of D/(M^3/(12c)) over M >= 8c:', float(worst))
