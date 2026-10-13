"""Print the non-binomial irreducible factors of the basis polynomials P_c (exact), with their exponent vectors
written as Laurent polynomials centred at the origin of the lattice (monomial factors stripped)."""
import sys, json
import flint
fn = sys.argv[1]; maxn = int(sys.argv[2])
D = json.load(open(fn)); A = D['A']; basis = D['basis']; d = len(A[0])
ctx = flint.fmpz_mpoly_ctx.get(('w', d), 'lex')
mins = [min(k[i] for k in A) for i in range(d)]
res = []
for c in sorted(basis, key=lambda c: sum(1 for x in c if x))[:maxn]:
    P = ctx.from_dict({tuple(k[i] - mins[i] for i in range(d)): x for k, x in zip(A, c) if x})
    for f, e in P.factor()[1]:
        if len(f.to_dict()) > 2:
            res.append(str(f))
for r in sorted(set(res)):
    print(r)
