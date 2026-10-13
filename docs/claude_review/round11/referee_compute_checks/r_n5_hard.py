# The two n=5 instances where the referee's zonotope search hit its limit: decide zeta = rho via the factorisation
# of an LLL basis of V(A, rho) (a basis polynomial that is monomial x two-term factors of total order rho is a
# zonotope certificate: its Newton polytope is a lattice zonotope inside conv(A)).
import math, itertools
import flint
from rzono import functionals, points, rank_poly
n = 5; base = 4; T = 3; L = 3
profs = [
 [[0, 3, 4], [0, 2, 3], [0, 3], [0, 2], [0, 4], [0, 2, 4], [0, 1, 3, 4], [0, 1, 3], [0]],
 [[0, 1, 4], [0, 2, 3], [0, 1, 2], [0, 2, 3, 4], [0, 1], [0, 2], [0, 4], [0, 2, 4], [0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 3]],
]
for prof in profs:
    prof = [frozenset(S) for S in prof]
    others, cons = functionals(n, prof, T, L, base)
    A = points(4, cons, T)
    m = 12
    mons = [al for al in itertools.product(range(2 * T + 1), repeat=4) if sum(al) < m]
    rows = [[math.prod(k[i] ** al[i] for i in range(4)) for k in A] for al in mons]
    M = flint.fmpz_mat(rows)
    r = M.rank()
    rows13 = [[math.prod(k[i] ** al[i] for i in range(4)) for k in A] for al in itertools.product(range(2*T+1), repeat=4) if sum(al) < m + 1]
    r13 = flint.fmpz_mat(rows13).rank()
    print(f"profile {[sorted(S) for S in prof]}: |A|={len(A)} rank Poly_<12={r} rank Poly_<13={r13} dimV(A,12)={len(A)-r}")
    N, nul = M.nullspace()
    B = flint.fmpz_mat([[int(N[i, j]) for i in range(len(A))] for j in range(nul)]).lll()
    ctx = flint.fmpz_mpoly_ctx.get(('w0', 'w1', 'w2', 'w3'), 'degrevlex')
    found = False
    for j in range(B.nrows()):
        c = [int(B[j, i]) for i in range(len(A))]
        P = ctx.from_dict({tuple(k[i] + T for i in range(4)): v for k, v in zip(A, c) if v})
        const, facs = P.factor()
        kinds = [(len(list(f.terms())), e, f.total_degree()) for f, e in facs]
        allbin = all(nt <= 2 for nt, e, dg in kinds)
        if allbin:
            ordr = sum(e for nt, e, dg in kinds if nt == 2)
            print(f"   basis poly {j}: product of monomial and two-term factors, total order {ordr}: ",
                  [str(f) + (f'^{e}' if e > 1 else '') for f, e in facs if len(list(f.terms())) == 2])
            found = found or ordr >= m
        else:
            print(f"   basis poly {j}: factor shapes (nterms, mult, deg): {kinds}")
    print("   zonotope certificate (zeta >= 12) found:", found)
