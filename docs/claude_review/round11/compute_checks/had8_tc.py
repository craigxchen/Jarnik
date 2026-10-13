"""Exact tangent-cone forms of the faithful had8 model and their structure.
Phi_c(u) = sum_k c_k <k,u>^m, c in V(A,m): exact integer coefficient vectors (core.phi_forms), LLL-reduced.
Reports: dim of the span, the set of variable-supports of monomials, minimal number of distinct variables per
monomial, and exact vanishing on candidate rational subspaces (coordinate planes, slab directions)."""
import sys, json, itertools
import flint
from core import phi_forms, monomials

fn = sys.argv[1]
D = json.load(open(fn))
A = [tuple(k) for k in D['A']]; basis = D['basis']; m = D['m']; d = len(A[0])
mons, rows = phi_forms(A, m, basis, d)
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
pr("forms: %d independent octics (degree %d) in %d variables; max |coef| after LLL = %d" % (
    len(rows), m, d, max(abs(x) for r in rows for x in r)))
supp = set()
minvars = 99
for r in rows:
    for a, c in zip(mons, r):
        if c:
            nv = sum(1 for x in a if x)
            minvars = min(minvars, nv)
            supp.add(tuple(i for i in range(d) if a[i]))
pr("min number of distinct variables in a monomial with nonzero coefficient:", minvars)
pr("number of distinct variable-supports:", len(supp), " sizes:", sorted(set(len(s) for s in supp)))
# coordinate subspaces span(e_S) contained in TC: all forms vanish when u_x = 0 for x not in S
def vanish_on_coord(S):
    for r in rows:
        for a, c in zip(mons, r):
            if c and all(a[i] == 0 for i in range(d) if i not in S):
                return False
    return True
for s in range(1, d + 1):
    good = [S for S in itertools.combinations(range(d), s) if vanish_on_coord(S)]
    pr("coordinate %d-planes in TC: %d of %d" % (s, len(good), len(list(itertools.combinations(range(d), s)))))
json.dump(dict(mons=mons, rows=rows, m=m, d=d), open(fn.replace("_basis.json", "_forms.json"), "w"))
open(fn.replace("_basis.json", "_tc_out.txt"), "w").write("\n".join(out) + "\n")
