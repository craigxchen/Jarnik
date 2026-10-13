"""Cone forms only for basis vectors whose non-binomial factors all have degree <= maxdeg (a subset suffices
for Lemma 5.1, since Bs is contained in the zero set of any subfamily)."""
import sys, json
import flint
import sympy as sp
fn = sys.argv[1]; maxdeg = int(sys.argv[2])
D = json.load(open(fn)); A = D['A']; basis = D['basis']; d = len(A[0])
ctx = flint.fmpz_mpoly_ctx.get(('w', d), 'lex')
mins = [min(k[i] for k in A) for i in range(d)]
X = sp.symbols('x0:%d' % d)
G = []
for c in basis:
    P = ctx.from_dict({tuple(k[i] - mins[i] for i in range(d)): x for k, x in zip(A, c) if x})
    g = sp.Integer(1); ok = True
    for f, e in P.factor()[1]:
        fd = f.to_dict()
        if len(fd) <= 2:
            if len(fd) == 2:
                assert sorted(int(v) for v in fd.values()) == [-1, 1]
            continue
        if f.total_degree() > maxdeg:
            ok = False; break
        fx = sp.expand(sum(int(cf) * sp.prod([(1 + X[i]) ** a[i] for i in range(d)]) for a, cf in fd.items()))
        assert len(set(sum(mm) for mm in sp.Poly(fx, *X).monoms())) == 1
        g = g * fx ** e
    if ok:
        G.append(sp.expand(g))
print("using", len(G), "of", len(basis), "basis vectors")
json.dump(dict(forms=[str(g) for g in G], d=d), open(fn.replace("_basis.json", "_sub_cone.json"), "w"))
