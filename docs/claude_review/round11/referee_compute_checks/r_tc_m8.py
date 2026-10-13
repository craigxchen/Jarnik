# Tangent-cone forms Phi_c(u) = sum_k c_k <k',u>^8 for the referee's own basis of V(A_2, 8) (M8), and the
# coordinate subspaces contained in TC (exact, flint).
import json, itertools
import flint
D = json.load(open('rk_had8_L2_m8.json'))
A, others, B, m = D['A'], D['others'], D['basis'], D['m']
U = flint.fmpz_mpoly_ctx.get(tuple(f"u{i}" for i in range(7)), 'degrevlex')
us = U.gens()
lin = [sum(k[x] * us[i] for i, x in enumerate(others)) for k in A]
pw = [l ** m for l in lin]
Phis = []
for c in B:
    P = U.from_dict({})
    for ck, p in zip(c, pw):
        if ck: P += ck * p
    Phis.append(P)
print("all Phi nonzero:", all(not P.is_zero() for P in Phis))
maxexp = max(max(e) for P in Phis for e in P.monoms())
minvars = min(sum(1 for t in e if t) for P in Phis for e in P.monoms())
print("max exponent of a variable in any Phi:", maxexp, "; min number of distinct variables in a monomial:", minvars)
print("max |coefficient|:", max(abs(int(c)) for P in Phis for c in P.coeffs()))
# coordinate subspaces span{e_i : i in S} contained in TC  <=>  every Phi vanishes when u_j = 0 for j not in S
# <=> no monomial of any Phi has support inside S.
supports = set(frozenset(i for i, t in enumerate(e) if t) for P in Phis for e in P.monoms())
for dim in range(1, 8):
    inTC = [S for S in itertools.combinations(range(7), dim) if not any(sp <= frozenset(S) for sp in supports)]
    print(f"dim {dim}: {len(inTC)} of {len(list(itertools.combinations(range(7), dim)))} coordinate subspaces lie in TC",
          (inTC if dim == 4 else ""))
# which 4-sets of non-base points are cuts (or complements of cuts) of H8 avoiding the base point 7?
from rmodels import sylvester_cuts
n, cuts = sylvester_cuts(3)
cut4 = sorted(tuple(sorted(S if 7 not in S else frozenset(range(8)) - S)) for S in cuts)
print("cuts (side not containing base 7):", cut4)
