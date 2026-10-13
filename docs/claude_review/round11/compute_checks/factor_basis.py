"""Factor the Laurent polynomials P_c(w) = sum_k c_k w^k (shifted to a polynomial) over Q for basis vectors c
of a model (exact, flint fmpz_mpoly).  Reports for each: number of terms, support size, factors (degree,
#terms, binomial or not), and whether each factor vanishes at e."""
import sys, json
import flint
fn = sys.argv[1]; maxn = int(sys.argv[2]) if len(sys.argv) > 2 else 10
D = json.load(open(fn)); A = D['A']; basis = D['basis']; d = len(A[0])
ctx = flint.fmpz_mpoly_ctx.get(('w', d), 'lex')
mins = [min(k[i] for k in A) for i in range(d)]
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
basis_sorted = sorted(basis, key=lambda c: sum(1 for x in c if x))
for c in basis_sorted[:maxn]:
    dct = {tuple(k[i] - mins[i] for i in range(d)): x for k, x in zip(A, c) if x}
    P = ctx.from_dict(dct)
    cont, facs = P.factor()
    desc = []
    for f, e in facs:
        fd = f.to_dict()
        if len(fd) == 1:
            continue  # monomial
        at_e = sum(int(v) for v in fd.values()) == 0
        desc.append("(deg %d, %d terms, %s, %s)^%d" % (f.total_degree(), len(fd), "binomial" if len(fd) == 2 else "NON-binomial",
                                                     "vanishes at e" if at_e else "nonzero at e", e))
        if len(fd) == 2:
            desc[-1] += " " + str(f)
    pr("support %d, sup %d, content %s: %s" % (len(dct), max(abs(x) for x in c), cont, "; ".join(desc)))
open(fn.replace("_basis.json", "_factor_out.txt"), "w").write("\n".join(out) + "\n")
