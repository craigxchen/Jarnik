"""T5: express the syl8 (L=2) e=6 invariant obstruction F (t2) through natural invariants, exactly (sympy).
Forms are reduced mod J(C) (drop monomials with an exponent >= 3) and written in the 8 orbit coordinates."""
import sys, pickle, itertools
import sympy as sp
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import syl
D = pickle.load(open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks/t2_syl8_L2_e6.pkl', 'rb'))
orbits = D['orbits']
k = sp.symbols('k0:7')
n, pats = syl(3)
sig = sum(k)
def lin(S):
    v = 0
    for x in S:
        v += (k[x] if x < 7 else -sig)
    return sp.expand(v)
s = [lin(S) for S in pats]
print("slab forms:", s)
def orbcoords(expr):
    P = sp.Poly(sp.expand(expr), *k)
    cd = {}
    for mon, c in P.terms():
        if max(mon) <= 2: cd[mon] = c
    out = []
    for orb in orbits:
        vals = [cd.get(a, 0) for a in orb]
        assert all(v == vals[0] for v in vals), "not invariant"
        out.append(vals[0])
    return out
q = sum(x**2 for x in k)
cands = {
  'sum s^6': sum(si**6 for si in s),
  'sum s^4 * sig^2': sum(si**4 for si in s) * sig**2,
  'sum s^2 * sig^4': sum(si**2 for si in s) * sig**4,
  'sig^6': sig**6,
  'sum s^4 * q': sum(si**4 for si in s) * q,
  'sum s^2 * q^2': sum(si**2 for si in s) * q**2,
  'q^3': q**3, 'q^2 sig^2': q**2 * sig**2, 'q sig^4': q * sig**4,
  'sum s^5 * sig': sum(si**5 for si in s) * sig,
  'sum s^3 * sig^3': sum(si**3 for si in s) * sig**3,
  'prod-type e6(k)': sum(sp.prod(c) for c in itertools.combinations(k, 6)),
  'sum_{p<q} s_p^3 s_q^3': sum(s[i]**3 * s[j]**3 for i in range(7) for j in range(i+1, 7)),
  'sum_{p<q} s_p^2 s_q^4+s_p^4 s_q^2': sum(s[i]**2 * s[j]**4 + s[i]**4*s[j]**2 for i in range(7) for j in range(i+1, 7)),
  'sum_{p<q<r} s_p^2 s_q^2 s_r^2': sum(s[i]**2 * s[j]**2 * s[t]**2 for i in range(7) for j in range(i+1, 7) for t in range(j+1, 7)),
  'sum s_p^4 s_p^... (s^2)^3 per p': 0,
}
F = [sp.Rational(x) for x in D['I'][0]]
JA = [[sp.Rational(x) for x in r] for r in D['VA']]
print("F (orbit coords):", F)
rows = {}
for name, ex in cands.items():
    if ex == 0: continue
    try:
        rows[name] = orbcoords(ex)
    except AssertionError:
        print("%-40s not invariant (sign)" % name); continue
    print("%-40s %s" % (name, rows[name]))
# which single candidate equals F mod J(A)^G (up to scalar)?
basisJ = sp.Matrix(JA)
for name, v in rows.items():
    M = sp.Matrix([list(F)] + JA + [v])
    Mj = sp.Matrix(JA + [v])
    if sp.Matrix(JA + [v]).rank() == 3 and M.rank() == 3:
        print("F == c *", name, "mod J(A)^G")
# span: F in span of candidates + J(A)?
names = list(rows)
Mall = sp.Matrix([rows[nm] for nm in names] + JA)
print("rank of candidates + J(A)^G:", Mall.rank(), " (8 = all invariants)")
# solve F = sum c_i cand_i + J part with few candidates: try pairs
import itertools as it
for r in (1, 2):
    for combo in it.combinations(names, r):
        M = sp.Matrix([rows[nm] for nm in combo] + JA)
        if sp.Matrix([rows[nm] for nm in combo] + JA + [F]).rank() == M.rank():
            sol = sp.Matrix([rows[nm] for nm in combo] + JA).T.solve_least_squares(sp.Matrix(F)) if False else None
            print("F in span of", combo, "+ J(A)^G")
