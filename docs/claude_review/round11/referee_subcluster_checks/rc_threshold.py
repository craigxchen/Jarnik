# Referee: where does the slab-box of [D] Thm 2.1 (T = n, w = n sqrt(8(n-1) ln 4r)) become LP-supercritical
# (m > 2w >= LP, hence V(A,m) != 0, dim >= |A|/2), and where do its integer generators become admissible in the
# sense of Lemma 1.1 (1.2) for W = W_min(n) = 3 n ln n - 14 n (two_thirds lower bound, C = 1)?
# Exact big-integer binomials for m; ln-bounds for |c|_1 as in [D] Thm 4.1: |c|_1 <= U^3 n^{2(m-1)}, U <= (2n+1)^{n-1}.
import math
def exact_m(n, T):
    rhs = (2*T+1)**(n-1)
    lo, hi = 0, 4*n*n*T
    while lo < hi:
        mid = (lo+hi+1)//2
        if 8*math.comb(mid+n-2, n-1) <= rhs: lo = mid
        else: hi = mid-1
    return lo
def stats(n, r):
    T = n; w = n*math.sqrt(8*(n-1)*math.log(4*r)); m = exact_m(n, T); D = 2*n*(n-1)
    lnc1 = 3*(n-1)*math.log(2*n+1) + 2*(m-1)*math.log(n)
    W = 3*n*math.log(n) - 14*n
    lhs = W*(m/2 - w)                   # need W(m/2 - LP/2) > 2 ln|c|_1 + 2m ln(D/2) - 2 ln m!, LP <= 2w
    rhs = 2*lnc1 + 2*m*math.log(D/2) - 2*math.lgamma(m+1)
    return m, w, (W > 0 and m > 2*w and lhs > rhs), lhs, rhs
first_sup = first_adm = first_4w = None
for n in range(64, 2400):
    m, w, adm, lhs, rhs = stats(n, n)
    if first_sup is None and m > 2*w: first_sup = (n, m, round(2*w, 1))
    if first_adm is None and adm: first_adm = (n, m, round(w, 1), round(lhs), round(rhs))
    if first_4w is None and m >= 4*w: first_4w = (n, m, round(4*w, 1))
print("r = n:")
print("  first n with m > 2w (slab-box LP-supercritical, V(A,m) != 0):", first_sup)
print("  first n with Lemma 1.1 admissibility of the generators at W = 3n ln n - 14n:", first_adm)
print("  first n with m >= 4w ([D]'s factor-2 margin):", first_4w)
print("  300 ln(4r) + 1 at r = n = 2137:", 300*math.log(4*2137)+1)
bad_sup = [n for n in range(440, 2400) if not exact_m(n, n) > 2*n*math.sqrt(8*(n-1)*math.log(4*n))]
bad_adm = [n for n in range(468, 2400) if not stats(n, n)[2]]
print("  n in [440,2400) with m <= 2w:", bad_sup[:10], "count", len(bad_sup))
print("  n in [468,2400) failing admissibility at W_min:", bad_adm[:10], "count", len(bad_adm))
for rr, lab in [(lambda n: int(n*math.log(n))+1, "r = n ln n"), (lambda n: n*n, "r = n^2")]:
    fs = next(n for n in range(64, 5000) if exact_m(n, n) > 2*n*math.sqrt(8*(n-1)*math.log(4*rr(n))))
    print(f"{lab}: first n with m > 2w: {fs}")
