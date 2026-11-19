"""Actual gcd (python-flint as gcd engine) of the five best5 raw marks at L=32, against fullgcd3."""
import time
from marks import raw_mark
from f9lib import F9P
from fams import fams
from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx
fp3 = fmpz_mod_poly_ctx(3); f9 = fq_default_ctx(modulus=fp3([1, 0, 1]), var='i'); R = fq_default_poly_ctx(f9)
ii = f9.gen()
forms, signs, g = fams['best5']; L = 32
t0 = time.time()
M = [raw_mark(l, e, 27, F9P([0]*L + [1]))[0] for l, e in zip(forms, signs)]
P = [R([int(a) + int(b)*ii for a, b in zip(z.re, z.im)]) for z in M]
G = P[0]
for x in P[1:]: G = G.gcd(x)
print('best5 L=32 actual gcd degree', G.degree(), '(fullgcd3 predicted 61470)', '%.0fs' % (time.time()-t0), flush=True)
Z = [p // G for p in P]
assert all(z * G == p for z, p in zip(Z, P))
D = Z[0].degree()
import itertools
dd = [(a - b).degree() for a, b in itertools.combinations(Z, 2)]
nm = [z * R([c**3 for c in z.coeffs()]) for z in Z]
print('D', D, 'diffs', dd, 'floor(D/2)', D//2, 'norms equal', all(x == nm[0] for x in nm), 'norm in F_3', all(c**3 == c for c in nm[0].coeffs()))
