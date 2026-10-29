"""T13: staircase (standard-monomial) sufficient criterion.  For a degree-compatible term order (deglex / degrevlex),
N(S) = standard monomials = greedy independent evaluation columns in increasing order (mod p).  A_l subset A gives
N(A_l) subset N(A).  If N(A)_e = Union_{l in L_pair} N(A_l)_e then (SHIFT_e) holds (dimension argument).
Reports, per degree e: |N(A)_e| = H_A(e), |Union_l N(A_l)_e|, and the uncovered monomials.  Usage: model L order"""
import sys, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
model, L, order = sys.argv[1], int(sys.argv[2]), sys.argv[3]
p = P1
n, pats = syl(3) if model == 'syl8' else sdc(7)
d = n - 1
A = build_box(n, [(S, L) for S in pats]); Aset = set(A)
mons = []
for e in range(0, 2 * d + 1):
    me = smonos(d, e)
    if order == 'deglex': me = sorted(me, reverse=True)                       # x0 > x1 > ... lex within degree
    elif order == 'degrevlex': me = sorted(me, key=lambda a: tuple(-v for v in reversed(a)), reverse=True)
    elif order == 'deglex_rev': me = sorted(me, key=lambda a: tuple(reversed(a)), reverse=True)
    mons += me
# increasing order = ascending degree, and within degree the order above taken from SMALLEST to largest:
# a term order: m' < m. We list within degree from smallest to largest by reversing the 'largest-first' lists.
order_list = []
for e in range(0, 2 * d + 1):
    me = [a for a in mons if sum(a) == e]
    order_list += list(reversed(me))
def staircase(S):
    E = evalmat(S, order_list)          # |mons| x |S|
    M = nmod(E.T.copy(), p)             # |S| x |mons|  columns = monomials in increasing order
    R, r = M.rref()
    piv = []
    for i in range(r):
        for j in range(R.ncols()):
            if int(R[i, j]) != 0: piv.append(j); break
    return set(order_list[j] for j in piv)
t0 = time.time()
NA = staircase(A)
print(model, "L=%d order=%s |A|=%d |N(A)|=%d (%.0fs)" % (L, order, len(A), len(NA), time.time() - t0), flush=True)
U = set()
for l in Lpair_half(d):
    U |= staircase(shiftset(A, l, Aset))
print("  (%.0fs)" % (time.time() - t0))
first_fail = None
for e in range(0, 2 * d + 1):
    Ne = [a for a in NA if sum(a) == e]
    if not Ne: break
    unc = [a for a in Ne if a not in U]
    print("  e=%d |N(A)_e|=%d uncovered=%d %s" % (e, len(Ne), len(unc), unc[:6]), flush=True)
    if unc and first_fail is None: first_fail = e
print("staircase criterion certifies SHIFT_e for e <", first_fail)
