"""T15: exact (over Q) confirmation for M8 (syl8, L = 2): sigma^6 notin J(A)_6, sigma^6 in J(A_l)_6 for l = e0, e0-e1
(by the range lemma and directly), sigma^7 in J(A)_7.  sigma = k_0 + ... + k_6 = -k_8."""
import sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
from math import factorial
n, pats = syl(3); d = 7
A = build_box(n, [(S, 2) for S in pats]); Aset = set(A)
def sigma_pow(e):
    F = {}
    for al in smonos(d, e):
        F[al] = factorial(e) // math.prod(factorial(a) for a in al)
    return F
for e in (6, 7):
    F = sigma_pow(e)
    print("e=%d: sigma^e in J(A)_e: %s" % (e, len(J_member_space(A, e, d, [F])) == 1))
F = sigma_pow(6)
for l in [(1, 0, 0, 0, 0, 0, 0), (1, -1, 0, 0, 0, 0, 0), (-1, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, -1)]:
    Al = shiftset(A, l, Aset)
    print("l=%s: sigma values on A_l: %s ; sigma^6 in J(A_l)_6: %s" % (l, sorted(set(sum(k) for k in Al)), len(J_member_space(Al, 6, d, [F])) == 1))
# rigidity of top level, exact enumeration
top = [k for k in A if sum(k) == 3]
print("top level sigma=3:", top)
