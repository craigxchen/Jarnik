"""Pell eight-point quartic only (see ref_curve_intersections.py for the method)."""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
sys.argv = [sys.argv[0]]
import ref_curve_intersections as R
from itertools import product
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit('/', 1)[0] + '/rerun')
from check_vojta_geometry import H_lin, F, pmul as hpmul, pconj as hpconj
from gauss import conj
words = []
for word in product((0, 1), repeat=4):
    p = sum(word)
    if p % 2 == 0: continue
    Pd = {(0, 0): F if p == 1 else conj(F)}
    for f, w in zip((H_lin(-2), H_lin(0), H_lin(2), H_lin(4)), word):
        Pd = hpmul(Pd, f if w else hpconj(f))
    P = [R.Z] * 5
    for (a, b), c in Pd.items():
        assert a + b == 4
        P[a] = R.gadd(P[a], (Fr(c[0]), Fr(c[1])))
    words.append(R.trim(P))
print("Pell eight-point quartic: 8 points, homogeneous degree 4 in (X:Y)")
R.report("Pell", words, 4, 7)
R.report("Pell", words, 4, 8)
diffs = [R.psub(p, words[0]) for p in words[1:]]
g = []
for f in diffs + [R.pconj(f) for f in diffs]:
    g = R.pgcd(g, f) if g else f
lc = R.ginv(g[-1]); g = [R.gmul(c, lc) for c in g]
print("   monic gcd of all 8-point differences (X, Y=1), coeffs low->high:", [(str(c[0]), str(c[1])) for c in g])
