"""T2: the syl8 (L=2) e=6 obstruction is a G-invariant: G = coordinate symmetry group of the box (perms of the 7
non-base coordinates preserving the slab set).  Exact over Q: invariant sextics (orbit sums of standard monomials) that
lie in J(A_l)_6 for l = e_0, e_0 - e_1 (orbit reps of L_pair; M_e(A_{-l}) = M_e(A_l)) but not in J(A)_6."""
import sys, time, itertools
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
L = int(sys.argv[1]) if len(sys.argv) > 1 else 2
e = int(sys.argv[2]) if len(sys.argv) > 2 else 6
n, pats = syl(3); d = 7
A = build_box(n, [(S, L) for S in pats]); Aset = set(A)
patset = set(frozenset(S) for S in pats) | set(frozenset(set(range(8)) - set(S)) for S in pats)
G = []
for perm in itertools.permutations(range(7)):
    pi = list(perm) + [7]
    if all(frozenset(pi[x] for x in S) in patset for S in pats):
        G.append(perm)
print("|G| =", len(G), flush=True)
# check box invariance
for g in G[:5]:
    assert all(tuple(k[g.index(x)] for x in range(7)) in Aset for k in A)
mons = smonos(d, e)
idx = {a: i for i, a in enumerate(mons)}
seen = set(); orbits = []
for a in mons:
    if a in seen: continue
    orb = set()
    for g in G:
        b = [0] * 7
        for x in range(7): b[g[x]] = a[x]
        orb.add(tuple(b))
    seen |= orb; orbits.append(sorted(orb))
print("standard monomials of degree", e, ":", len(mons), " G-orbits:", len(orbits), flush=True)
forms = [{a: 1 for a in orb} for orb in orbits]
t0 = time.time()
VA = J_member_space(A, e, d, forms)
print("dim (J(A)_e)^G =", len(VA), "(%.0fs)" % (time.time() - t0), flush=True)
reps = {'e0': (1, 0, 0, 0, 0, 0, 0), 'e0-e1': (1, -1, 0, 0, 0, 0, 0)}
spaces = {}
for name, l in reps.items():
    Al = shiftset(A, l, Aset)
    spaces[name] = J_member_space(Al, e, d, forms)
    print("dim (J(A_%s)_e)^G =" % name, len(spaces[name]), " |A_l| =", len(Al), "(%.0fs)" % (time.time() - t0), flush=True)
# intersection of the two invariant spaces (exact, via fmpq nullspace)
def intersect(U, V):
    import sympy
    if not U or not V: return []
    m = len(U[0])
    Uq = [[sympy.Rational(int(x.p), int(x.q)) for x in r] for r in U]
    Vq = [[sympy.Rational(int(x.p), int(x.q)) for x in r] for r in V]
    M = sympy.Matrix([[Uq[i][j] for i in range(len(U))] + [-Vq[i][j] for i in range(len(V))] for j in range(m)])
    out = []
    for c in M.nullspace():
        out.append([sum(c[i] * Uq[i][j] for i in range(len(U))) for j in range(m)])
    if not out: return []
    R = sympy.Matrix(out).rref()[0]
    return [[flint.fmpq(int(sympy.fraction(R[i, j])[0]), int(sympy.fraction(R[i, j])[1])) for j in range(m)] for i in range(R.rank())]
I = intersect(spaces['e0'], spaces['e0-e1'])
print("dim (cap_l J(A_l)_e)^G =", len(I), "  vs dim (J(A)_e)^G =", len(VA))
def rank(vs):
    return 0 if not vs else flint.fmpq_mat(vs).rank()
print("J(A)^G subset cap:", rank(I + VA) == len(I))
# a representative of the obstruction: vector in I not in span VA; reduce against VA
basis = list(VA)
obs = [v for v in I if rank(basis + [v]) > len(basis)]
if obs:
    v = obs[0]
    print("obstruction representative (orbit coefficients):")
    for c, orb in zip(v, orbits):
        if c != 0: print("   coef", c, " orbit size", len(orb), " rep", orb[0])
import pickle
pickle.dump({'orbits': orbits, 'VA': [[str(x) for x in r] for r in VA], 'I': [[str(x) for x in r] for r in I], 'G': G},
            open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks/t2_syl8_L%d_e%d.pkl' % (L, e), 'wb'))
