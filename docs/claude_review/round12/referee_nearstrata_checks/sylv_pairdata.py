# Concrete Sylvester profile (k=16, n=65536, columns 1[a.x=0], a!=0, e_p=1) with genuine primes p_a = 1 mod 4
# in a short interval [P0, P0+H].  N = prod p_a (C <= 1).  Compute exactly-ish (float64 + error bound) the pair
# distances d_xy = d(x xor y) via a Walsh-Hadamard transform, the pair excesses eps = d - W/2, abar, and test
# (i) Gamma(binary point)=2.9996 <= abar  (so averaged (MIN) fails at window orders by Prop 5.2), and
# (ii) the individual pair gaps |u_x-u_y| >= exp(-eps_xy/2) are infeasible for n points in [-1,1]
#      (sum of n-1 sorted adjacent gaps would exceed 2), so the profile-specific pair data exclude every u.
import numpy as np, math
k = 16; n = 2**k
P0 = int(math.exp(25)); H = 6_000_000
lo, hi = P0, P0 + H
lim = int(math.isqrt(hi)) + 1
small = np.ones(lim+1, bool); small[:2] = False
for i in range(2, int(lim**0.5)+1):
    if small[i]: small[i*i::i] = False
sp = np.nonzero(small)[0]
seg = np.ones(hi-lo+1, bool)
for p in sp:
    st = max(p*p, ((lo + p - 1)//p)*p)
    seg[st-lo::p] = False
primes = lo + np.nonzero(seg)[0]
primes = primes[primes % 4 == 1]
print('primes = 1 mod 4 in [P0, P0+H]:', len(primes), ' P0 =', P0)
assert len(primes) >= n-1
pr = primes[:n-1]
logs = np.log(pr.astype(np.float64))
print('log p range:', logs.min(), logs.max(), 'spread', logs.max()-logs.min())
# weights indexed by a in F_2^k \ 0
wvec = np.zeros(n); wvec[1:] = logs
W = wvec.sum()
# Walsh-Hadamard: S(z) = sum_a (-1)^{a.z} w_a ; d(z) = (W - S(z))/2
S = wvec.copy(); h = 1
while h < n:
    S = S.reshape(-1, 2*h)
    a = S[:, :h].copy(); b = S[:, h:].copy()
    S[:, :h] = a + b; S[:, h:] = a - b
    S = S.reshape(-1); h *= 2
d = (W - S)/2
eps = d[1:] - W/2
abar = W/(2*(n-1))
print('W =', W, ' 3n ln n - 13.91 n =', 3*n*math.log(n) - 13.91*n, ' T3 consistent:', W >= 3*n*math.log(n) - 13.91*n)
print('abar =', abar, ' (>= Gamma(binary point) = 2.99962):', abar >= 2.99962)
print('eps_xy: min', eps.min(), 'max', eps.max(), ' (Lemma 3.1 needs eps>=0:', eps.min() >= 0, ')')
thr = 2*math.log((n-1)/2)
print('max eps =', eps.max(), ' vs 2 ln((n-1)/2) =', thr, ' -> individual gaps infeasible:', eps.max() < thr - 1e-6)
print('min allowed gap exp(-max eps/2) =', math.exp(-eps.max()/2), ' times (n-1) =', (n-1)*math.exp(-eps.max()/2), '(>2 means infeasible)')
