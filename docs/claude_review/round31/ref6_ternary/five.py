"""Exact t-degrees of the specialised raw differences (F_i - F_j)(t, t^L) for q=27, k=3 families (own F_9[t]
arithmetic with FFT multiplication). Certified cluster: z_j = F_j/G(t,t^L), D = (n-g)L,
deg(z_i - z_j) = deg(F_i - F_j) - gL <= D/2 required."""
import sys, itertools, time
from marks import raw_mark
from f9lib import F9P
from fams import fams
name = sys.argv[1]; Ls = [int(x) for x in sys.argv[2:]]
forms, signs, g = fams[name]
q = 27; k = 3; n = q**k; h = n//q
for L in Ls:
    t0 = time.time()
    T = F9P([0]*L + [1])
    M = [raw_mark(l, e, q, T)[0] for l, e in zip(forms, signs)]
    assert all(z.deg() == n*L and z.lead() == (1, 0) for z in M)
    D = (n-g)*L
    dd = [(M[a]-M[b]).deg() - g*L for a, b in itertools.combinations(range(len(M)), 2)]
    assert min(dd) >= 0   # distinct
    ok = 2*max(dd) <= D
    print('%s L=%d D=(n-g)L=%d parent=%d  max reduced diff=%d  D/2=%.1f  cluster(with formal gcd only): %s  [%.0fs]'
          % (name, L, D, 2*D, max(dd), D/2, ok, time.time()-t0), flush=True)
    print('   reduced diff degrees', dd, flush=True)
