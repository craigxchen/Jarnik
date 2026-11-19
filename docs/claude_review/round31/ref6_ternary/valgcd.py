"""Validate fullgcd3 against an actual gcd computation at L=3 (python-flint used only as gcd engine here)."""
from marks import raw_mark
from f9lib import F9P
from fams import fams
from fullgcd3 import extra3
from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx
fp3 = fmpz_mod_poly_ctx(3); f9 = fq_default_ctx(modulus=fp3([1, 0, 1]), var='i'); R = fq_default_poly_ctx(f9)
ii = f9.gen()
for name in ['astra5', 'best5']:
  for L in [3, 4]:
    forms, signs, g = fams[name]
    T = F9P([0]*L + [1])
    M = [raw_mark(l, e, 27, T)[0] for l, e in zip(forms, signs)]
    P = [R([int(a) + int(b)*ii for a, b in zip(z.re, z.im)]) for z in M]
    G = P[0]
    for x in P[1:]: G = G.gcd(x)
    res = extra3(forms, signs, L)
    print(name, L, 'actual gcd deg', G.degree(), ' predicted g*L+extra', g*L + list(res.values())[0][1], flush=True)
