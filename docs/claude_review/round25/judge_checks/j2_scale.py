# Judge check J2: (a) exact rational check of referee Prop R1(d); (b) where the R1 and cube data sit
# relative to the growth regime W = O(M log M); (c) referee coset finding 1: unit upper-triangular example.
from fractions import Fraction as F
from math import log, comb
print("(a) R1(d): (1-2^-m)(1+eta)/(1-eta) < (1-eta) with m=2t+1, eta=2^-(2t+4)")
for t in range(3, 31):
    m = 2*t+1; eta = F(1, 2**(2*t+4))
    lhs = (1 - F(1, 2**m))*(1+eta)/(1-eta); rhs = 1-eta
    assert lhs < rhs, t
print("   holds exactly for t = 3..30")
print("(a') route-5 regime n-1 >= 300 ln(4r), r = 2^(2t+1)-2, M = largest 4q+1 <= 2^t")
for t in range(10, 16):
    r = 2**(2*t+1)-2; M = 2**t - 3  # 2^t - 3 = 1 mod 4 for t>=2
    print(f"   t={t}: M={M}, M-1={M-1}, 300 ln 4r = {300*log(4*r):.1f}, ok={M-1 >= 300*log(4*r)}")
print("(b) scale of the data (lower bounds on W = log N'):")
t = 13; m = 2*t+1; r = 2**m - 2; M = 8189
w0 = (4*t+5)*log(2); eta = 2.0**(-(2*t+4))
W_R1 = r*(1-eta)*w0
print(f"   R1 data t=13: r={r}, w0>={w0:.2f}, W>={W_R1:.3e}; M log M = {M*log(M):.3e}; "
      f"ratio W/(M log M) >= {W_R1/(M*log(M)):.3e}")
print(f"   generally r = 2^(2t+1)-2 >= 2M^2-2 primes, so W >= (2M^2-2) log 5; e.g. M=8189: {(2*M*M-2)*log(5):.3e}")
for Mc in (10, 20, 40):
    Wc = (2**(Mc-1)-1)*log(5)
    print(f"   cube class M={Mc}: W >= (2^(M-1)-1) log 5 = {Wc:.3e}; M log M = {Mc*log(Mc):.1f}")
print("(c) U = I - (strict upper ones), r x r: exact inverse")
def inv_unit_upper(r):
    U = [[F(1) if i == j else (F(-1) if j > i else F(0)) for j in range(r)] for i in range(r)]
    # back substitution for U X = I
    X = [[F(0)]*r for _ in range(r)]
    for col in range(r):
        for i in reversed(range(r)):
            s = F(1 if i == col else 0) - sum(U[i][k]*X[k][col] for k in range(i+1, r))
            X[i][col] = s
    return X
for r in (5, 10, 20):
    X = inv_unit_upper(r)
    mx = max(abs(x) for row in X for x in row)
    rowsum = max(sum(abs(x) for x in row) for row in X)
    ok = all(X[i][j] == (2**(j-i-1) if j > i else (1 if i == j else 0)) for i in range(r) for j in range(r))
    print(f"   r={r}: entries = 2^(j-i-1) above diagonal: {ok}; max entry {mx}; max row l1 {rowsum} = 2^(r-1)")
print("(d) Lemma R2 constant: 1 + log m <= (1.39 - (log 4 - 1)) m for m = 1..10^6:",
      all(1+log(k) <= (1.39-(log(4)-1))*k for k in range(1, 10**6+1)))
