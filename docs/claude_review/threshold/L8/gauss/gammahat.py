"""Circle height class gamma-hat on Mbar_{0,k+2} (labels P,Q, circle points 1..k):
gamma.D_S = 1/4 if {S,S^c} = {{P,Q},[k]};  2^-k if S = {P} u A (A nonempty proper subset of [k]);  0 otherwise.
(1) verify all Keel relations exactly; (2) k=4: test against all 15 Keel-Vermeire divisors (Hassett-Tschinkel:
Eff(Mbar_{0,6}) = cone(boundary, KV)) under every labelling."""
from fractions import Fraction as Fr
from itertools import combinations, permutations

def splits(labels):
    labels = list(labels); n = len(labels); first = labels[0]; out = []
    for r in range(2, n - 1):
        for S in combinations(labels, r):
            S = frozenset(S); Sc = frozenset(labels) - S
            if first in S and len(Sc) >= 2:
                out.append(S)
            elif first not in S and len(S) >= 2 and r < n - r + 0 and False:
                pass
    # canonical representative: the side containing labels[0]
    reps = set()
    for r in range(2, n - 1):
        for S in combinations(labels, r):
            S = frozenset(S); Sc = frozenset(labels) - S
            if len(Sc) >= 2:
                reps.add(S if first in S else Sc)
    return sorted(reps, key=lambda s: (len(s), sorted(map(str, s))))

def gamma(k):
    labels = ['P', 'Q'] + list(range(1, k + 1)); circ = frozenset(range(1, k + 1))
    g = {}
    for S in splits(labels):
        Sc = frozenset(labels) - S
        if S == frozenset(['P', 'Q']) or Sc == frozenset(['P', 'Q']):
            g[S] = Fr(1, 4)
        elif ('P' in S and 'Q' not in S) or ('Q' in S and 'P' not in S):
            g[S] = Fr(1, 2 ** k)
        else:
            g[S] = Fr(0)
    return labels, g

def keel_ok(labels, g):
    labels = list(labels); bad = 0
    for i, j, k_, l in permutations(labels, 4):
        if not (i < j if all(isinstance(x, int) for x in (i, j)) else str(i) < str(j)):
            continue
        def side(S, a):
            return S if a in S else frozenset(labels) - S
        lhs = sum(v for S, v in g.items() if side(S, i) >= {i, j} and not (side(S, i) & {k_, l}))
        rhs = sum(v for S, v in g.items() if side(S, i) >= {i, k_} and not (side(S, i) & {j, l}))
        if lhs != rhs:
            bad += 1
    return bad == 0

for k in range(2, 8):
    labels, g = gamma(k)
    print('k=%d: Keel relations hold for gamma-hat: %s' % (k, keel_ok(labels, g)))

# ---- k = 4: Keel-Vermeire divisors in Kapranov psi_6 model: KV = 2H - sum_{i<=5} E_i - E_13 - E_14 - E_23 - E_24
# for partition {12|34|56} with 6 the psi-point.  Dictionary: E_i = D_{i6}, E_ij = D_{ij6},
# H = D_{lm} + E_i + E_j + E_k + E_ij + E_ik + E_jk  for {i,j,k,l,m} = {1..5}.
def KV_dot(gd):  # gd: function frozenset(S) of {1..6} -> gamma.D_S (S any side)
    E = lambda *a: gd(frozenset(a) | {6})
    H = gd(frozenset({4, 5})) + E(1) + E(2) + E(3) + E(1, 2) + E(1, 3) + E(2, 3)
    return 2 * H - sum(E(i) for i in range(1, 6)) - E(1, 3) - E(1, 4) - E(2, 3) - E(2, 4)

labels, g = gamma(4)
def make_gd(assign):  # assign: dict 1..6 -> label in labels
    def gd(S):
        T = frozenset(assign[x] for x in S); Tc = frozenset(labels) - T
        return g[T] if T in g else g[Tc]
    return gd
vals = set()
for perm in permutations(labels):
    assign = dict(zip(range(1, 7), perm))
    vals.add(KV_dot(make_gd(assign)))
print('k=4: gamma-hat . KV over all labellings:', sorted(vals))
# sanity: the KV formula must give the same number for a curve class symmetric... check with beta (all 1):
print('sanity beta.KV (beta.D_S=1 for all S):', KV_dot(lambda S: Fr(1)))
