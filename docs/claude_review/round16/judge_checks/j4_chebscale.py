# Judge check J4 (50-digit arithmetic): the referee's 1-D Chebyshev witness refutes the cube form AC_T(.;1/2)
# for every u unless W > W_cheb := 4 m' ln(sqrt(d)/0.275) + 4 ln(2(m'+1)/(C d)) (C=1, |u|_2 <= sqrt d worst case).
# Compare with W_ref = 3 n ln n (near-extremal scale of the record) at regime sizes, r = n.
import mpmath as mp
mp.mp.dps=50
for n in (2966, 4096, 10**4, 10**5, 10**6, 10**8):
    d=n-1; r=n
    w=mp.sqrt(8*d*mp.log(4*r)); m=int(mp.ceil(4*w))
    regime = d >= 316*mp.log(4*r)
    Wc=4*m*mp.log(mp.sqrt(d)/mp.mpf('0.275'))+4*mp.log(2*(m+1)/mp.mpf(d))
    Wr=3*n*mp.log(n)
    print("n=%d regime=%s m'=%d W_cheb=%.4g 3nln n=%.4g ratio=%.3f"%(n,regime,m,float(Wc),float(Wr),float(Wc/Wr)))
# Part 2: smallest admissible W for actual clusters (two_thirds, C=1): W_lo = 3 n ln n - 14 n.
# Find the n-range (regime, r=n) where W_lo < W_cheb, i.e. the cube form is false for EVERY u at the
# smallest admissible W.  Note: median 0.55|u|_2 is the Gaussian approximation used by the referee (evidence-level
# constant); the bisection below is on that formula.
def Wc(n):
    d=n-1; r=n; w=mp.sqrt(8*d*mp.log(4*r)); m=int(mp.ceil(4*w))
    return 4*m*mp.log(mp.sqrt(d)/mp.mpf('0.275'))+4*mp.log(2*(m+1)/mp.mpf(d))
last=None
for n in range(2966, 40000):
    if 3*n*mp.log(n)-14*n < Wc(n): last=n
print("largest n in [2966,40000) with 3n ln n - 14n < W_cheb:", last)
for n in (2966, 5000, 8000, last, last+1):
    print(n, float(3*n*mp.log(n)-14*n), float(Wc(n)))
