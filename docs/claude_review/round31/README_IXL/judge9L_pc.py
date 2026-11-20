# Judge check (iii): the printed lattice form (L_8^pc) of expert_note.md (lines 216-222) and Bal^pc.
# Own code, exact integers.
import itertools, random, math
random.seed(31)
out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)
pairs = list(itertools.combinations(range(7), 2))
def v_p(n, p):
    if n == 0: return 99
    v = 0
    while n % p == 0: n //= p; v += 1
    return v
# Lattice form, local at p: pairwise coprime e_T (3<=|T|<=6) => at most one e_T divisible by p.
# c_ij divisible by p iff v_p(y_i - y_j) > v_p(prod_{T' >= {i,j}} e_T') = v*[T >= {i,j}].
def min_divisible(p, v, Tsize):
    T = set(range(Tsize)) if Tsize else set()
    best = 99
    mod = p**(v+1)
    if Tsize:
        # T-marks congruent mod p^v: WLOG class 0, lifts j*p^v (j<p); outside marks arbitrary mod p^{v+1}
        Tchoices = [k*p**v for k in range(p)]
        rng = [Tchoices if i in T else range(mod) for i in range(7)]
    else:
        rng = [range(p)]*7
    for y in itertools.product(*rng):
        c = 0
        for (i, j) in pairs:
            den = v if (Tsize and i in T and j in T) else 0
            if v_p(y[i]-y[j], p) > den if (y[i]-y[j]) % mod else True:
                # difference divisible by mod (> den always since den <= v < v+1) or valuation > den
                c += 1
        best = min(best, c)
        if best <= 1: break
    return best
for p in (2, 3):
    say(f"p={p}, no block divisible by p: min #cofactors divisible by p = {min_divisible(p, 0, 0)}")
    for Tsize in (3, 4, 5, 6):
        for v in ((1, 2, 3) if p == 2 else (1,)):
            if p == 3 and v == 2 and Tsize == 3:
                pass
            say(f"p={p}, one block |T|={Tsize} with v_p(e_T)={v}: min #cofactors divisible by p = {min_divisible(p, v, Tsize)}")
# Bal^pc: visible T (2<=|T|<=6) active at p (a residue class of x mod p^t, t>=1)
def active(x, p, tmax=40):
    act = set()
    for t in range(1, tmax):
        cls = {}
        for i, xi in enumerate(x): cls.setdefault(xi % p**t, set()).add(i)
        for c in cls.values():
            if 2 <= len(c) <= 6: act.add(frozenset(c))
        if all(len(c) == 1 for c in cls.values()): break
    return act
from collections import Counter
cnt2, cnt3 = Counter(), Counter(); big2 = 0
for _ in range(20000):
    while True:
        x = [0] + random.sample(range(-10**6, 10**6), 6)
        g = 0
        for xi in x: g = math.gcd(g, xi - x[0])
        if g == 1 and len(set(x)) == 7: break
    a2 = active(x, 2); a3 = active(x, 3)
    cnt2[len(a2)] += 1; cnt3[len(a3)] += 1
    if sum(1 for T in a2 if len(T) >= 3) < 2: big2 += 1
say("Bal^pc: #visible T active at 2 over 20000 primitive 7-tuples:", dict(cnt2), "; tuples with <2 active sets of size>=3:", big2)
say("Bal^pc: #visible T active at 3:", dict(sorted(cnt3.items())))
# corrected statement is locally satisfiable: referee example
y = (0, 5, 10, 107, 267, 69, -229)
e = {frozenset({0, 1, 2}): 5}
def strip23(n):
    n = abs(n)
    for p in (2, 3):
        while n % p == 0: n //= p
    return n
cof = {}
okdiv = True
for (i, j) in pairs:
    den = math.prod(v for T, v in e.items() if i in T and j in T)
    if (y[j]-y[i]) % den: okdiv = False
    cof[(i, j)] = (y[j]-y[i])//den
cop = all(math.gcd(strip23(cof[a]), strip23(cof[b])) == 1 for a, b in itertools.combinations(pairs, 2))
say("corrected form example y=(0,5,10,107,267,69,-229), e_{123}=5: divisibility", okdiv, "; cofactors pairwise coprime outside {2,3}:", cop)
open(__file__.replace('.py', '.out'), 'w').write("\n".join(out) + "\n")
