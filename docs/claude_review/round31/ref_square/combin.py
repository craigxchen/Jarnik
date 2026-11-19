"""Referee's own exact/combinatorial checks for Astra sections 3-4."""
import random, itertools, math
from itertools import combinations, product
from fractions import Fraction
from scipy.optimize import linprog
import numpy as np

random.seed(7)
marks = range(5)
pairs = list(combinations(marks, 2))
# 15 cuts as sides not containing mark 4 (index 4 = "mark five")
cuts = [frozenset(s) for r in range(1, 5) for s in combinations(range(4), r)]
assert len(cuts) == 15


def cross(A, e):
    return (e[0] in A) != (e[1] in A)


# --- each pair separated by 8 cuts
for e in pairs:
    assert sum(cross(A, e) for A in cuts) == 8

# --- (P) pairing average lemma: exhaustive on [0,e]^4, e<=7
for e in range(1, 8):
    for a in product(range(e+1), repeat=4):
        s = abs(a[0]+a[1]-a[2]-a[3]) + abs(a[0]+a[2]-a[1]-a[3]) + abs(a[0]+a[3]-a[1]-a[2])
        assert s <= 3*e
print("pairing lemma: sum of three matchings <= 3e_p on [0,e]^4, e<=7: OK")

# --- (5): |t^2-1| >= sqrt2 D(t^2)^(-1/4) for random t = A/conj(A) (exact-ish via high precision)
from mpmath import mp, mpf, mpc, fabs, sqrt
mp.dps = 60
ps = [5, 13, 17, 29, 37, 41, 53, 61]
pis = {5: (2, 1), 13: (3, 2), 17: (4, 1), 29: (5, 2), 37: (6, 1), 41: (5, 4), 53: (7, 2), 61: (6, 5)}
worst = mpf(10)
for trial in range(20000):
    A = mpc(1, 0)
    D = 1
    for p in ps:
        k = random.randint(-3, 3)
        if k == 0:
            continue
        x, y = pis[p]
        g = mpc(x, y)/mpc(x, -y)
        A *= g**k
        D *= p**abs(k)
    eps = [mpc(1, 0), mpc(0, 1), mpc(-1, 0), mpc(0, -1)][random.randint(0, 3)]
    t = eps*A
    s = t*t
    if fabs(s-1) < mpf(10)**-50:
        continue
    r = fabs(s-1)/(sqrt(2)*mpf(D*D)**(-0.25))
    worst = min(worst, r)
    assert r >= 1-mpf(10)**-40
print("(5) on 20000 random squares: min ratio", mp.nstr(worst, 8))

# --- common-mean lemma: P,Q disjoint, sum P x = sum Q x  => E|X-Y| <= 1/2 ; adversarial search
best = 0
for trial in range(20000):
    n1 = random.randint(1, 3); n2 = random.randint(1, 4-n1+1)
    P = np.random.dirichlet(np.ones(n1)); Q = np.random.dirichlet(np.ones(n2))
    # maximise sum P_i Q_j |x_i - y_j| subject to P.x = Q.y, x,y in [0,1] -- vertices: random sign pattern LP
    for _ in range(3):
        sg = np.random.choice([-1, 1], size=(n1, n2))
        # linear objective given ordering signs: sum PQ sg (x_i - y_j); maximise under constraints
        c = -np.concatenate([ (P[:, None]*Q[None, :]*sg).sum(1), -(P[:, None]*Q[None, :]*sg).sum(0)])
        res = linprog(c, A_eq=[np.concatenate([P, -Q])], b_eq=[0], bounds=[(0, 1)]*(n1+n2))
        if res.status == 0:
            x = res.x[:n1]; y = res.x[n1:]
            val = sum(P[i]*Q[j]*abs(x[i]-y[j]) for i in range(n1) for j in range(n2))
            best = max(best, val)
assert best <= 0.5 + 1e-9
print("common-mean lemma: max E|X-Y| found", round(best, 6), "(bound 1/2)")

# --- cut decomposition inequality and costs 12, 8
for trial in range(5000):
    e = random.randint(1, 6)
    a = [random.randint(0, e) for _ in marks]
    a[random.randint(0, 4)] = 0
    a[random.choice([i for i in marks if a[i] == 0] and list(marks))] = e
    if min(a) != 0 or max(a) != e:
        continue
    c = [random.randint(-3, 3) for _ in range(4)]
    c.append(-sum(c))
    lhs = abs(sum(ci*ai for ci, ai in zip(c, a)))
    rhs = 0
    for t in range(1, e+1):
        A = frozenset(i for i in marks if a[i] >= t)
        rhs += abs(sum(c[i] for i in A))
    assert lhs <= rhs
for c, val in [((1, 1, -1, -1, 0), 12), ((1, -1, 0, 0, 0), 8)]:
    assert sum(abs(sum(c[i] for i in A)) for A in cuts) == val
    # also with mark 5 in the word: all orientations equal by sum zero
for perm in itertools.permutations(range(5)):
    c = [0]*5
    c[perm[0]] = c[perm[1]] = 1; c[perm[2]] = c[perm[3]] = -1
    assert sum(abs(sum(c[i] for i in A)) for A in cuts) == 12
print("cut decomposition inequality, costs 12 (all four-mark words) and 8: OK")


# --- (7) and the squarefree claims via LP over pair distributions
def lp_value(S):
    """min over distributions mu on the 10 pairs of max_{A in S} mu(cross A)."""
    nP = len(pairs)
    c = np.zeros(nP+1); c[-1] = 1
    A_ub = []; b_ub = []
    for A in S:
        row = [1.0 if cross(A, e) else 0.0 for e in pairs] + [-1.0]
        A_ub.append(row); b_ub.append(0)
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=[[1.0]*nP+[0]], b_eq=[1], bounds=[(0, None)]*(nP+1))
    return res.fun


two = [A for A in cuts if len(A) == 2 or len(A) == 3 and False]
# two-mark cuts: |A| = 2 or |A| = 3 (complement = {x,4})
twomark = [A for A in cuts if len(A) in (2, 3)]
single = [A for A in cuts if len(A) in (1, 4)]
assert len(twomark) == 10 and len(single) == 5
# (7) general: for each two-mark cut E, uniform on 7 pairs meeting E's 2-side
for E in twomark:
    side = E if len(E) == 2 else frozenset(range(5)) - E
    Qs = [e for e in pairs if e[0] in side or e[1] in side]
    assert len(Qs) == 7
    for A in cuts:
        k = sum(cross(A, e) for e in Qs)
        assert (k == 6) if A == E else (k <= 4)
print("(7) certificate: 7 pairs, cut E crosses 6, others <= 4: OK")
v_all_but_one = max(lp_value([A for A in cuts if A != E]) for E in twomark)
print("LP value with one two-mark cut missing (all other 14 present):", round(v_all_but_one, 6), " 4/7 =", round(4/7, 6))
# at most six cuts
worst6 = 0
for S in combinations(cuts, 6):
    worst6 = max(worst6, lp_value(S))
print("max LP value over all 6-cut supports:", round(worst6, 6), "(claim <= 1/2)")
worst7 = 0; arg7 = None
for S in combinations(cuts, 7):
    v = lp_value(S)
    if v > worst7:
        worst7, arg7 = v, S
print("max LP value over all 7-cut supports:", round(worst7, 6), [sorted(A) for A in arg7])
worst9 = 0
for S in combinations(cuts, 9):
    worst9 = max(worst9, lp_value(S))
print("max LP value over all 9-cut supports:", round(worst9, 6), "(claim <= 4/7 =", round(4/7, 6), ")")
print("LP value all 15 cuts:", round(lp_value(cuts), 6), " (3/5 = Cilleruelo-Cordoba)")
print("LP value ten two-mark cuts:", round(lp_value(twomark), 6))

# --- angle model: lattice, M(c)
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
B0 = Matrix([[int(i in e) - int(4 in e) for e in pairs] for i in range(4)])
print("Smith form of anchored 4x10 pair matrix:", [smith_normal_form(B0)[i, i] for i in range(4)])


def M(c):
    return sum(abs(c[i]+c[j]) for i, j in pairs)


cnt = {}
for c4 in product(range(-4, 5), repeat=4):
    c = list(c4) + [-sum(c4)]
    if all(x == 0 for x in c):
        continue
    m = M(c)
    n1 = sum(abs(x) for x in c)
    assert 2*m >= 3*n1
    isroot = n1 == 2
    if isroot:
        assert m == 6
    else:
        assert m >= 8, c
    cnt[m] = cnt.get(m, 0) + 1
print("M(c) over c in [-4,4]^4: roots have M=6, all others M>=8; smallest values:", sorted(cnt.items())[:5])
