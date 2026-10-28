# Referee J4, check 5 (exact): explicit constants of the proof.
#   eta(delta) := (1/(1+2delta) - 15/16)/2  (> 0 iff delta < 1/30); kappa := 1 - (15/16+eta)(1+2delta) > 0
#   K := max(ceil(9+2/eta), M_1),  M_1 := least M0 with (eta/2) C(M,5) > D M^4 for all M >= M0 (checked on a range)
#   R_1(C) := max(C^(2/kappa), (C R_0)^2)   [Lemma 1.2: R'' >= R/chord >= sqrt(R)/C >= R_0]
# Also: general k, the convexity-averaging threshold (1 - 2^{1-k})(1+2delta) < 1  <=>  delta < 1/(2^k-2).
from fractions import Fraction as Fr
from math import comb, ceil
for delta in (Fr(1, 31), Fr(1, 40), Fr(1, 64)):
    eta = (1/(1+2*delta) - Fr(15, 16))/2
    kappa = 1 - (Fr(15, 16)+eta)*(1+2*delta)
    for D in (1, 2, 5, 20):
        # g(M) = C(M,5)/M^4 = M prod_{i=1..4}(1-i/M)/120 is increasing in M >= 5, so bisect for the first M with
        # (eta/2) C(M,5) > D M^4; verify exactly at M1-1 and on [M1, M1+2000].
        ok = lambda M: eta/2*comb(M, 5) > D*M**4
        lo, hi = 9, 10**12
        while lo < hi:
            mid = (lo+hi)//2
            if ok(mid): hi = mid
            else: lo = mid+1
        M1 = lo
        assert (M1 == 9 or not ok(M1-1)) and all(ok(M) for M in range(M1, M1+2000))
        K = max(ceil(9 + 2/eta), M1)
        print("delta=%s eta=%s kappa=%s D=%d: M_1=%d (<= 937.3 D/eta = %.0f), K=%d" %
              (delta, eta, kappa, D, M1, 937.3*D/eta, K))
for k in range(3, 9):
    thr = Fr(1, 2**k - 2)
    lhs_at = (1 - Fr(2)**(1-k))*(1+2*thr)
    print("k=%d: (1-2^{1-k})(1+2/(2^k-2)) = %s  (=1 exactly: threshold is 1/(2^k-2))" % (k, lhs_at))
