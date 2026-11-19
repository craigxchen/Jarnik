"""Minimal L with the FULL gcd (exact, via fullgcd3) for q=27, k=3 families; raw differences by own FFT arithmetic."""
import sys, itertools, time
from marks import raw_mark
from f9lib import F9P
from fams import fams
from fullgcd3 import extra3
name = sys.argv[1]; Ls = [int(x) for x in sys.argv[2:]]
forms, signs, g = fams[name]
q = 27; n = q**3
for L in Ls:
    t0 = time.time()
    T = F9P([0]*L + [1])
    M = [raw_mark(l, e, q, T)[0] for l, e in zip(forms, signs)]
    raw = [(M[a]-M[b]).deg() for a, b in itertools.combinations(range(len(M)), 2)]
    res = extra3(forms, signs, L)
    g0, ex = list(res.values())[0]; assert g0 == g
    gam = g*L + ex
    D = n*L - gam
    md = max(raw) - gam
    print('%s L=%d formal gL=%d extra=%d full gcd=%d  D=%d parent=%d  max diff=%d  floor(D/2)=%d  cluster=%s [%.0fs]'
          % (name, L, g*L, ex, gam, D, 2*D, md, D//2, md <= D//2, time.time()-t0), flush=True)
