"""Exact check of Lemma 2.2 (theory.md): for a primitive integral 7-configuration x (anchor at infinity),
the congruence lattice
   Lambda(x) = { y in Z^7 : y_i = y_j mod p^t whenever {i,j} lie in a residue class C of x mod p^t, |C| >= 3 }
has index  prod_{3<=|T|<=6} d_T^{|T|-1}  (d_T = contact integers), and
   log|x_i - x_j| = sum_{T visible, T ∋ i,j} log d_T.
Contacts computed from residue classes (independent of the index computation, which uses a Smith normal form)."""
import random, itertools, math
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form
from sympy import factorint

rng = random.Random(7)
def planted(m=7):
    # x_i = sum over a few primes of (cluster-structured multiples)
    x = [rng.randint(0, 10**3) for _ in range(m)]
    for p in rng.sample([3, 5, 7, 11, 13], 3):
        e = rng.randint(1, 3)
        S = rng.sample(range(m), rng.randint(3, 6))
        base = rng.randint(0, 10**3)
        for i in S: x[i] = base + (p**e) * rng.randint(0, 50)
    return x

def contacts(x):
    m = len(x); ps = set()
    for i, j in itertools.combinations(range(m), 2): ps |= set(factorint(abs(x[i]-x[j])))
    d = {}; cls = {}
    for p in ps:
        t = 1
        while True:
            q = p**t; groups = {}
            for i, v in enumerate(x): groups.setdefault(v % q, []).append(i)
            big = [frozenset(g) for g in groups.values() if len(g) >= 2]
            if not big: break
            for g in big:
                if len(g) < m: d[g] = d.get(g, 1) * p
                if len(g) >= 3: cls.setdefault((p, t), []).append(g)
            t += 1
    return d, cls

ok = 0; tested = 0
while tested < 60:
    x = planted()
    if len(set(x)) < 7: continue
    g = 0
    for i in range(7): g = math.gcd(g, x[i] - x[0])
    if g != 1: continue
    tested += 1
    d, cls = contacts(x)
    # identity log|x_i - x_j| = sum d_T
    ok1 = all(abs(x[i]-x[j]) == math.prod(v for T, v in d.items() if i in T and j in T)
              for i, j in itertools.combinations(range(7), 2))
    pred = math.prod(v**(len(T)-1) for T, v in d.items() if len(T) >= 3)
    rows = []; mods = []
    for (p, t), gs in cls.items():
        for C in gs:
            C = sorted(C)
            for i in C[1:]:
                r = [0]*7; r[i] = 1; r[C[0]] = -1; rows.append(r); mods.append(p**t)
    if rows:
        A = Matrix(rows); r = len(rows)
        Mfull = A.row_join(Matrix.diag(*mods))
        S = smith_normal_form(Mfull, domain=ZZ)
        diag = [abs(S[i, i]) for i in range(r)]
        idx = math.prod(mods) // math.prod(diag)
    else:
        idx = 1
    ok += (ok1 and idx == pred)
    if not (ok1 and idx == pred): print("MISMATCH", x, idx, pred, ok1)
print(f"tested {tested} primitive planted configurations: {ok} pass (identity (1.1) and index formula)")
