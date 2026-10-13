"""Exact cone reduction of ZE for the had8 models.
For each basis vector c_j of V(A,m): factor P_j over Q into a monomial, binomials (w^a - w^b) and non-binomial
irreducible factors f.  Check exactly that every non-binomial f satisfies ord_e f = deg f, i.e. that
f(1 + x) is a FORM (homogeneous polynomial) in x = w - 1.  Then
   Bs(A,m) minus the union of the subtori {w^{a-b} = 1}  is contained in  { w : x = w - 1 in Cone(G_1..G_s) },
G_j(x) = prod of the non-binomial factors of P_j, written in x (a form of degree m - ord_e(binomial part)).
Outputs the forms G_j (as sympy expressions in x0..x{d-1}) to <model>_cone.json."""
import sys, json
import flint
import sympy as sp
fn = sys.argv[1]
D = json.load(open(fn)); A = D['A']; basis = D['basis']; d = len(A[0]); m = D['m']
ctx = flint.fmpz_mpoly_ctx.get(('w', d), 'lex')
mins = [min(k[i] for k in A) for i in range(d)]
X = sp.symbols('x0:%d' % d)
G = []
Gform = []
allhom = True
summary = {}
binoms = set()
for c in basis:
    P = ctx.from_dict({tuple(k[i] - mins[i] for i in range(d)): x for k, x in zip(A, c) if x})
    cont, facs = P.factor()
    g = sp.Integer(1)
    ordb = 0; shapes = []
    for f, e in facs:
        fd = f.to_dict()
        if len(fd) == 1:
            continue
        if len(fd) == 2:
            vals = sorted(int(v) for v in fd.values())
            assert vals == [-1, 1], "two-term factor not of the form w^a - w^b"
            ordb += e
            binoms.add(str(f))
            continue
        fx = sp.expand(sum(int(cf) * sp.prod([(1 + X[i]) ** a[i] for i in range(d)]) for a, cf in fd.items()))
        poly = sp.Poly(fx, *X)
        degs = set(sum(mon) for mon in poly.monoms())
        homog = (len(degs) == 1)
        allhom = allhom and homog
        shapes.append((f.total_degree(), len(fd), e, min(degs), "form" if homog else "NOT-form"))
        g = g * fx ** e
    key = (ordb, tuple(sorted(shapes)))
    summary[key] = summary.get(key, 0) + 1
    G.append(g)
    if all(sh[4] == "form" for sh in shapes):
        Gform.append(g)
    assert ordb + sum(s[3] * s[2] for s in shapes) == m
print("model", fn, "m =", m, "dimV =", len(basis))
print("shape of each basis polynomial: (order of binomial part, sorted (deg, #terms, mult) of non-binomial factors): count")
for k, v in summary.items():
    print("   ", k, v)
print("all non-binomial factors are forms in x = w - 1:", allhom, "; order of binomial part + orders of the other factors = m (checked)")
print("distinct binomial factors:", len(binoms), sorted(binoms)[:12], "...")
print("basis vectors whose non-binomial part is a form in x:", len(Gform), "of", len(G))
print("all binomial factors:", sorted(binoms))
json.dump(dict(forms=[str(sp.expand(g)) for g in Gform], d=d, homogeneous=True, n_used=len(Gform), n_basis=len(G)), open(fn.replace("_basis.json", "_cone.json"), "w"))
