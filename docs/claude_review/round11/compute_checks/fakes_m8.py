"""Random (fake-type) positions delta ~ U[0,Delta]^n for a model: certify sum_k c_k e^{i<k,delta>} != 0 for some
basis vector c, by the leading Taylor term beating the remainder (50 digits):
   |sum_k c_k <k,u>^m/m!| > sum_k |c_k| |<k,u>|^{m+1}/(m+1)! * e^{|<k,u>|}   (u_x = delta_x - delta_n)
since e^{ix} - sum_{t<=m}(ix)^t/t! is bounded by |x|^{m+1}/(m+1)! e^{|x|} and the terms t < m vanish (c in V)."""
import sys, json, math
import numpy as np
import mpmath as mp
mp.mp.dps = 50
fn = sys.argv[1]; trials = int(sys.argv[2]); Delta = mp.mpf(sys.argv[3]); seed = int(sys.argv[4])
D = json.load(open(fn)); A = D['A']; basis = D['basis']; m = D['m']; d = len(A[0])
rng = np.random.default_rng(seed)
cert = 0
for t in range(trials):
    delta = [mp.mpf(float(x)) * Delta for x in rng.uniform(0, 1, size=d + 1)]
    u = [delta[x] - delta[d] for x in range(d)]
    xs = [mp.fsum(k[i] * u[i] for i in range(d)) for k in A]
    ok = False
    for c in basis:
        lead = abs(mp.fsum(ck * x ** m for ck, x in zip(c, xs))) / mp.factorial(m)
        rem = mp.fsum(abs(ck) * abs(x) ** (m + 1) * mp.exp(abs(x)) for ck, x in zip(c, xs) if ck) / mp.factorial(m + 1)
        if lead > rem:
            ok = True
            break
    cert += ok
s = "%s: %d/%d random positions delta ~ U[0,%s]^%d certified OFF Z_P (leading Taylor term > remainder bound)" % (
    fn, cert, trials, mp.nstr(Delta, 3), d + 1)
print(s)
open(fn.replace("_basis.json", "_fakes_out.txt"), "w").write(s + "\n")
