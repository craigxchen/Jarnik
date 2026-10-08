# Referee: exact minimum of D = sum_{tau=1}^{e} (S_tau - M/2)^2 over level vectors n_0..n_e with
# n_0>=1, n_e>=1, 0<=n_a<=c, sum = M, any e>=1 (S_tau = n_0+...+n_{tau-1}).
# 4D = sum (2S-M)^2 is an integer.  S_1=n_0 in [1,c]; steps in [0,c]; last S_e in [M-c, M-1].
# Zero steps only repeat a term (cost >= 0), so WLOG steps in [1,c]: shortest path DP.
from fractions import Fraction

def minD4(M, c):
    INF = 10**30
    best = [INF] * (M + 1)
    for s in range(1, min(c, M - 1) + 1):
        best[s] = (2 * s - M) ** 2
    for s in range(1, M):
        if best[s] == INF:
            continue
        for t in range(s + 1, min(s + c, M - 1) + 1):
            v = best[s] + (2 * t - M) ** 2
            if v < best[t]:
                best[t] = v
    return min(best[s] for s in range(max(1, M - c), M))

if __name__ == '__main__':
    viol = 0; worst = None; tested = 0; worst_any = None
    for c in range(1, 41):
        for M in range(4 * c + 1, 4 * c + 400):
            d4 = minD4(M, c); tested += 1
            bound = Fraction((M - 4 * c) ** 3, 12 * c)
            r_any = Fraction(d4, 4) / bound
            worst_any = r_any if worst_any is None or r_any < worst_any else worst_any
            if Fraction(d4, 4) < bound:
                viol += 1; print('VIOLATION', M, c, d4 / 4, float(bound))
            if M >= 8 * c:
                r = Fraction(d4, 4) / Fraction(M ** 3, 12 * c)
                worst = r if worst is None or r < worst else worst
    print('tested', tested, 'violations', viol)
    print('min D/((M-4c)^3/(12c)) over M>4c:', float(worst_any))
    print('min D/(M^3/(12c)) over M>=8c:', float(worst))
    for c, M in ((1, 3000), (4, 3000), (12, 3000), (40, 4000)):
        print('c', c, 'M', M, 'Dmin/(M^3/(12c)) =', minD4(M, c) / 4 / (M ** 3 / (12 * c)))
