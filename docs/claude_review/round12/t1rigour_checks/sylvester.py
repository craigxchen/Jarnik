"""(A) Recheck of J1 (round11/README.md §5.1): column identity, exact m at T = n, cheapness inequality.
(B) Sylvester-Hadamard profile at T = 1: when are the k = log2 n coordinate columns cheap together
    (2 j k < m'), i.e. when does the braid statement fail at T = 1 (t1rigour.md §6.3)."""
import math, itertools
import mpmath
mpmath.mp.dps = 60

# (A1) column identity: for a != b nonzero in F_2^k and each class c of (a.x, b.x), find integers
# alpha, beta, gamma, delta (delta = +-2) with 1[a.x=0] = alpha + beta 1[b.x=0] + gamma 1[(a+b).x=0] + delta 1[(a.x,b.x)=c].
def dot(a, x):
    return bin(a & x).count('1') & 1
fails = 0; checked = 0
for k in range(3, 8):
    n = 1 << k
    # solve once per class pattern: values on the four classes determine the identity
    sols = {}
    for c in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        found = None
        for al, be, ga, de in itertools.product(range(-3, 4), repeat=4):
            if abs(de) != 2:
                continue
            good = True
            for ax, bx in [(0, 0), (0, 1), (1, 0), (1, 1)]:
                lhs = 1 if ax == 0 else 0
                rhs = al + be * (1 if bx == 0 else 0) + ga * (1 if (ax ^ bx) == 0 else 0) + de * (1 if (ax, bx) == c else 0)
                if lhs != rhs:
                    good = False; break
            if good:
                found = (al, be, ga, de); break
        sols[c] = found
    for a in range(1, n):
        for b in range(1, n):
            if a == b:
                continue
            for c, (al, be, ga, de) in sols.items():
                for x in range(n):
                    checked += 1
                    lhs = 1 if dot(a, x) == 0 else 0
                    rhs = al + be * (dot(b, x) == 0) + ga * (dot(a ^ b, x) == 0) + de * ((dot(a, x), dot(b, x)) == c)
                    if lhs != rhs:
                        fails += 1
print("J1 column identity: classes ->", sols, f"; checked {checked} (a,b,c,x) for k=3..7, failures {fails}")

def m_Tn(n):
    """largest m with 8 binom(m+n-2, n-1) <= (2n+1)^(n-1)  (T = n)."""
    d = n - 1
    rhs = d * mpmath.log(2 * n + 1) - mpmath.log(8)
    lo, hi = 1, 3 * n * n
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mpmath.loggamma(mid + n - 1) - mpmath.loggamma(n) - mpmath.loggamma(mid) <= rhs:
            lo = mid
        else:
            hi = mid - 1
    return lo

for k in (13, 14):
    n = 1 << k
    m = m_Tn(n)
    exact = None
    if k == 13:
        R = (2 * n + 1) ** (n - 1)
        ok1 = 8 * math.comb(m + n - 2, n - 1) <= R
        ok2 = 8 * math.comb(m + 1 + n - 2, n - 1) <= R
        exact = (ok1, not ok2)
    w = n * math.sqrt(8 * (n - 1) * math.log(4 * (n - 1)))
    cost = 4 * math.floor(w / 2) + (n // 4) * 2 * n
    print(f"J1 T=n, n=2^{k}: m = {m} (exact boundary check: {exact}), 0.7n(n-1) = {0.7*n*(n-1):.6g}, "
          f"cost 4floor(w/2) + n^2/2 = {cost:.6g}, cheap (cost < m): {cost < m}, cost < 0.7n(n-1): {cost < 0.7*n*(n-1)}")
# monotonicity of the ratio (n/2 + 2 sqrt(8(n-1) ln 4(n-1)))/(0.7 (n-1)) for n >= 2^14 (upper bound for cost/(0.7n(n-1)))
prev = None
for k in range(14, 41):
    n = 2.0 ** k
    rat = (n / 2 + 2 * math.sqrt(8 * (n - 1) * math.log(4 * (n - 1)))) / (0.7 * (n - 1))
    assert prev is None or rat < prev
    prev = rat
print("J1 ratio bound decreasing on 2^14..2^40 and < 1 at 2^14:",
      (2**14 / 2 + 2 * math.sqrt(8 * (2**14 - 1) * math.log(4 * (2**14 - 1)))) / (0.7 * (2**14 - 1)))

def m_T1(n):
    d = n - 1
    rhs = d * mpmath.log(3) - mpmath.log(8)
    lo, hi = 1, 2 * n
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mpmath.loggamma(mid + n - 1) - mpmath.loggamma(n) - mpmath.loggamma(mid) <= rhs:
            lo = mid
        else:
            hi = mid - 1
    m = lo
    if n <= 70000:
        assert 8 * math.comb(m + n - 2, n - 1) <= 3 ** (n - 1) < 8 * math.comb(m + n - 1, n - 1)
    return m

print("\nSylvester at T = 1 (r = n - 1): k, n, w, j, 2jk, ceil(4w), m_max, m_max>=4w, k cols cheap at m_max (2jk < m_max)")
for k in range(11, 23):
    n = 1 << k
    d = n - 1
    w = math.sqrt(8 * d * math.log(4 * (n - 1)))
    j = math.floor(w / 2)
    mm = m_T1(n)
    print(f"{k}, {n}, {w:.1f}, {j}, {2*j*k}, {math.ceil(4*w)}, {mm}, {mm >= 4*w}, {2*j*k < mm}")
