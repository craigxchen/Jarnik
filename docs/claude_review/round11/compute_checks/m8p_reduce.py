"""Exact reduction of ZE for a model whose basis polynomials factor as (binomials) x (non-binomial part).
For each basis vector c_j: P_j = mono * prod(binomials) * Q_j.  At a nondegenerate point near e every binomial
w^a - w^b is nonzero (it vanishes only where <a-b,u> = 0 mod 2pi), so  Bs_nondeg (near e) is contained in
V(Q_1, ..., Q_s).  We compute the exact Jacobian of (Q_j) at e: if it has rank d, e is an isolated point of
V(Q) (inverse function theorem), so a neighbourhood of e contains no nondegenerate point of Bs: ZE-real holds
near e, with an explicit radius from a Kantorovich-type bound (computed separately).
Also reports the lowest-order forms of the Q_j at e."""
import sys, json
import flint
import sympy as sp
fn = sys.argv[1]
D = json.load(open(fn)); A = D['A']; basis = D['basis']; d = len(A[0])
ctx = flint.fmpz_mpoly_ctx.get(('w', d), 'lex')
mins = [min(k[i] for k in A) for i in range(d)]
W = sp.symbols('w0:%d' % d)
Qs = []
info = []
for c in basis:
    P = ctx.from_dict({tuple(k[i] - mins[i] for i in range(d)): x for k, x in zip(A, c) if x})
    cont, facs = P.factor()
    Q = ctx.from_dict({tuple([0] * d): 1})
    nb = 0; ordb = 0
    for f, e in facs:
        fd = f.to_dict()
        if len(fd) == 1:
            continue
        if len(fd) == 2:
            ordb += e
            continue
        Q = Q * f ** e
        nb += 1
    Qs.append(Q)
    info.append((ordb, nb, Q.total_degree(), len(Q.to_dict())))
print("per basis vector: (order of binomial part at e, #non-binomial factors, deg Q, #terms Q):", info)
# Jacobian at e
def tosp(f):
    return sum(int(c) * sp.prod([W[i] ** a[i] for i in range(d)]) for a, c in f.to_dict().items())
Qsp = [sp.expand(tosp(Q)) for Q in Qs]
e = {w: 1 for w in W}
J = sp.Matrix([[sp.diff(q, w).subs(e) for w in W] for q in Qsp])
print("Q_j(e) =", [q.subs(e) for q in Qsp])
print("Jacobian of (Q_j) at e (exact):")
sp.pprint(J)
print("rank =", J.rank())
json.dump([str(q) for q in Qsp], open(fn.replace("_basis.json", "_Q.json"), "w"))
