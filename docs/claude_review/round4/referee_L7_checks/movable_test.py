"""EXACT (referee): necessary conditions for the class cub13 of L7 Prop 3 to be the class of a COVERING
family of curves.  A covering family has nonnegative degree on every effective divisor.  Tests:
 (1) pushforward to every Mbar_{0,6} (forget one label) against all 15 Keel-Vermeire divisors (pulled back
     KV divisors are effective on Mbar_{0,7});
 (2) pushforward to every Mbar_{0,5} (forget two labels) must have self-intersection >= 0 (nef on a surface).
KV formula (L8 6.3, psi_6 model, pairing (12)(34)(56)):
   KV = 2H - sum_{i<=5}E_i - E_13 - E_14 - E_23 - E_24,  H = D_45+E_1+E_2+E_3+E_12+E_13+E_23,
   E_i = D_{i6}, E_ij = D_{ij6}.  Checked: beta_6.KV = 5, forgetful fibre .KV in {0,1}; independence of the
   label assignment within the pairing is checked on Keel-valid classes.
"""
import sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0, '../L7')
import classes7 as C

def key(S, labels):
    """canonical key for split S of 'labels': frozenset side not containing max(labels)."""
    S = frozenset(S); top = max(labels)
    return S if top not in S else frozenset(labels) - S

def push(vec, labels, forget):
    """vec: dict on splits of 'labels' (keys by key()); returns dict on splits of labels-forget."""
    new = [x for x in labels if x != forget]
    out = {}
    for r in range(2, len(new) - 1):
        for S in itertools.combinations(new, r):
            k = key(S, new)
            if k in out: continue
            a = vec[key(S, labels)]; b = vec[key(set(S) | {forget}, labels)]
            out[k] = a + b
    return out, new

def KV(vec6, labels, assign):
    """assign: tuple (l1..l6) mapping model labels 1..6 -> actual labels; pairing (l1l2)(l3l4)(l5l6)."""
    L = dict(zip(range(1, 7), assign))
    D = lambda *xs: vec6[key([L[x] for x in xs], labels)]
    E = lambda *xs: D(*(xs + (6,)))
    H = D(4, 5) + E(1) + E(2) + E(3) + E(1, 2) + E(1, 3) + E(2, 3)
    return 2*H - sum(E(i) for i in range(1, 6)) - E(1, 3) - E(1, 4) - E(2, 3) - E(2, 4)

def pairings(labels):
    labels = list(labels)
    if not labels: yield []; return
    a = labels[0]
    for b in labels[1:]:
        rest = [x for x in labels if x not in (a, b)]
        for pr in pairings(rest): yield [(a, b)] + pr

def kv_values(vec6, labels):
    vals = {}
    for pr in pairings(labels):
        vs = set()
        for perm in itertools.permutations(pr):           # order of pairs
            for flips in itertools.product([0, 1], repeat=3):
                assign = []
                for (x, y), f in zip(perm, flips): assign += ([x, y] if not f else [y, x])
                vs.add(KV(vec6, labels, tuple(assign)))
        vals[tuple(pr)] = vs
    return vals

# sanity on Mbar_{0,6}: beta_6 and a forgetful fibre
lab6 = [1, 2, 3, 4, 5, 6]
beta6 = {key(S, lab6): 1 for r in (2, 3) for S in itertools.combinations(lab6, r)}
v = kv_values(beta6, lab6)
print('beta_6 . KV values:', set().union(*v.values()))
beta7 = {S: 1 for S in C.SPL}
f7 = {S: (1 if (C.is_pair(S) and len(S) == 5) else 0) for S in C.SPL}    # fibre of forgetting 7? (D_{a7})
for j in [7, 1]:
    p6, new = push(f7, list(range(1, 8)), j)
    print('fibre f_7 pushed (forget %d): KV values' % j, set().union(*kv_values(p6, new).values()))

planes = list(itertools.combinations(range(1, 7), 3))
best = None
for choice in itertools.combinations(planes, 13):
    vec = C.kapranov(3, {}, {}, {p: 1 for p in choice})
    if min(vec.values()) >= 0: best = choice; break
cub = C.kapranov(3, {}, {}, {p: 1 for p in best})

lab7 = list(range(1, 8))
neg = []
for j in lab7:
    p6, new = push(cub, lab7, j)
    vals = kv_values(p6, new)
    for pr, vs in vals.items():
        assert len(vs) == 1, (j, pr, vs)
        x = vs.pop()
        if x < 0: neg.append((j, pr, x))
print('cub13: (pi_j)_* cub13 . KV < 0 for', len(neg), 'of', 7*15, '(j, pairing) pairs')
for t in neg[:10]: print('   forget', t[0], 'pairing', t[1], 'value', t[2])
# surface test
bad = []
for j, k in itertools.combinations(lab7, 2):
    p6, new = push(cub, lab7, j); p5, new5 = push(p6, new, k)
    top = max(new5)
    m = {a: p5[key([a, top], new5)] for a in new5 if a != top}
    # e = psi_top degree on Mbar_{0,5}: psi_n = sum_{S containing n, j0,k0 not in S} D_S
    others = [a for a in new5 if a != top]; j0, k0 = others[:2]
    e = sum(p5[key(S, new5)] for r in (2, 3) for S in itertools.combinations(new5, r)
            if top in S and j0 not in S and k0 not in S and len(S) <= 3)
    sq = e*e - sum(x*x for x in m.values())
    bad.append(((j, k), e, sq))
print('Mbar_{0,5} projections (forget j,k): min self-intersection', min(b[2] for b in bad),
      ' values:', sorted(set(b[2] for b in bad)))
