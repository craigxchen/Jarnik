# Judge check J4b (rigorous constants, 50-digit arithmetic).
# Claim: for every u with |u|_inf = 1 (so |u|_2^2 <= d), X = <k',u>, k' uniform on {-1,0,1}^d:
#   Var X = (2/3)|u|_2^2, so Pr(|X| < s) >= 1/2 for s = sqrt(4d/3)  (Chebyshev's inequality).
#   Q = 2^(1-m) s^m T_m(X/s) has leading form X^m and |Q| <= 2^(1-m) s^m on {|X| < s}.
# Hence AC_T(u;1/2) fails for EVERY such u as soon as 2^(1-m) s^m <= d^m eps_0, eps_0 = C d e^(-W/4)/(m+1), i.e.
#   W <= W_rig := 4 m ln(sqrt(3d)) + 4 ln(C d/(2(m+1)))      (C = 1 below).
# Compare with W_lo = 3 n ln n - 14 n (two_thirds, C = 1), r = n.
import mpmath as mp
mp.mp.dps=50
def par(n):
    d=n-1; r=n; w=mp.sqrt(8*d*mp.log(4*r)); m=int(mp.ceil(4*w)); return d,m
def Wrig(n):
    d,m=par(n); return 4*m*mp.log(mp.sqrt(3*d))+4*mp.log(mp.mpf(d)/(2*(m+1)))
# sanity: exact check of the Chebyshev bound on [-1,1] for small m
for m in (3,7,12):
    T=mp.chebyt
    mx=max(abs(T(m,mp.mpf(x)/1000)) for x in range(-1000,1001))
    print("max|T_%d| on grid = %s"%(m,mp.nstr(mx,10)))
last=None
for n in range(2966, 20000):
    d,m=par(n)
    if 3*n*mp.log(n)-14*n <= Wrig(n): last=n
print("largest n in [2966,20000) with 3n ln n - 14n <= W_rig:", last)
for n in (2966, 3500, last, last+1, 10**4, 10**5):
    d,m=par(n); print("n=%d m'=%d W_lo=%.6g W_rig=%.6g 3nlnn=%.6g"%(n,m,float(3*n*mp.log(n)-14*n),float(Wrig(n)),float(3*n*mp.log(n))))
