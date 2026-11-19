# Exact dichotomy computations: given five integer columns V_i, integers n_i, blocks, compute everything.
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd, prod, log
import sympy
edges = list(combinations(range(5), 2))
CYC = {(0,1),(1,2),(2,3),(3,4),(0,4)}
def det(u, v): return u[0]*v[1]-u[1]*v[0]
def lcm(a, b): return a*b//gcd(a, b)

def analyse(V, J, F, label, verbose=True):
    """V: five integer 2-vectors; J: dict edge->block; F: complemented edges."""
    S = [{e for e in edges if (i in e) != (e in F)} for i in range(5)]
    n = [prod(J[e] for e in S[i]) for i in range(5)]
    for e, f in combinations(edges, 2): assert gcd(J[e], J[f]) == 1, "blocks not private"
    g = {(i,j): prod(J[e] for e in S[i] & S[j]) for i, j in edges}
    for i, j in edges: assert gcd(n[i], n[j]) == g[i,j]
    m = {}; t = {}
    for i, j in edges:
        m[i,j] = det(V[i], V[j]); m[j,i] = -m[i,j]
        assert m[i,j] % g[i,j] == 0, "divisibility fails"
        t[i,j] = m[i,j]//g[i,j]; t[j,i] = -t[i,j]
        assert t[i,j] != 0
    prim = all(gcd(*v) == 1 for v in V)
    delta = [V[i][0]**2+V[i][1]**2-n[i] for i in range(5)]
    Ht = prod(abs(t[e]) for e in edges)
    h = 4 if not F else 5
    # (IV.5) left sides, exactly
    c = [Fr(1, prod(t[i,j] for j in range(5) if j != i)) for i in range(5)]
    lhs = [sum(c[i]*t[i,k]*J[tuple(sorted((i,k)))] for i in range(5) if i != k) for k in range(5)]
    rhs = [-sum(c[i]*t[i,k]*J[tuple(sorted((i,k)))]*Fr(delta[i], n[i]) for i in range(5) if i != k) for k in range(5)]
    assert lhs == [-r for r in [-x for x in rhs]] or True
    for k in range(5): assert lhs[k] == rhs[k]
    den_ok = all((lhs[k]*prod(abs(t[a,b]) for a,b in edges if k not in (a,b))).denominator == 1 for k in range(5))
    # Q_*: solve on three columns, check on all five
    Mq = sympy.Matrix([[V[i][0]**2, V[i][0]*V[i][1], V[i][1]**2] for i in range(3)])
    sol = Mq.LUsolve(sympy.Matrix(n[:3]))
    qs = [sympy.Rational(x) for x in sol]
    Qv = [qs[0]*v[0]**2+qs[1]*v[0]*v[1]+qs[2]*v[1]**2 for v in V]
    interp_ok = all(Qv[i] == n[i] for i in range(5))
    den = sympy.ilcm(*[x.q for x in qs])
    E = [1-qs[0], -qs[1], 1-qs[2]]
    detE = E[0]*E[2]-E[1]**2/4
    # Gram check for every pair
    Ebil = lambda u, w: E[0]*u[0]*w[0]+E[1]*(u[0]*w[1]+u[1]*w[0])/2+E[2]*u[1]*w[1]
    for i, j in edges:
        if not interp_ok: break
        assert Ebil(V[i], V[i]) == delta[i]
        assert delta[i]*delta[j]-Ebil(V[i], V[j])**2 == m[i,j]**2*detE
    md = max(abs(d) for d in delta)
    Elist = [abs(Ebil(V[i],V[j])) for i,j in edges]
    out = dict(label=label, h=h, prim=prim, delta=delta, signs=sorted({(d>0)-(d<0) for d in delta}),
               t=[t[e] for e in edges], Ht=Ht, lhs_zero=all(x == 0 for x in lhs), lhs_den_ok=den_ok,
               interp_ok=interp_ok, den=den, den_divides_Ht=(Ht % den == 0), detE=detE,
               detE_lower_ok=(detE == 0 or abs(detE) >= sympy.Rational(1, 4*Ht**2)),
               rank=(0 if all(x == 0 for x in E) else (1 if detE == 0 else 2)), maxdelta=md, maxE=max(Elist))
    if out['rank'] == 1:
        # find lambda, primitive l with E = lambda l^2
        M = sympy.Matrix([[E[0], E[1]/2], [E[1]/2, E[2]]])
        col = M[:,0] if M[:,0] != sympy.zeros(2,1) else M[:,1]
        l = [sympy.Rational(x) for x in col]; L = sympy.ilcm(*[x.q for x in l]); l = [x*L for x in l]
        G = sympy.igcd(*[int(x) for x in l]); l = [int(x)//G for x in l]
        lam = (E[0]/l[0]**2) if l[0] != 0 else (E[2]/l[1]**2)
        assert all(delta[i] == lam*(l[0]*V[i][0]+l[1]*V[i][1])**2 for i in range(5))
        out['lambda'] = lam; out['ell'] = l
        out['lam_ge_1/Ht'] = abs(lam) >= sympy.Rational(1, Ht)
    if verbose:
        P = (prod(J.values()))**(1/10) if False else None
        logP = log(prod(J.values()))/10
        kap = log(max(abs(x) for x in out['t']))/logP
        out['logP_maxdelta'] = (log(md)/logP) if md else None
        out['kappa'] = kap
        out['threshold_exp'] = h-3-12*kap
        out['detE_abs_log_P'] = (log(abs(float(detE)))/logP) if detE else None
        for k in ['label','h','prim','signs','lhs_zero','lhs_den_ok','interp_ok','den','den_divides_Ht','rank','detE_lower_ok',
                  'kappa','threshold_exp','logP_maxdelta','detE_abs_log_P','lambda','ell','lam_ge_1/Ht']:
            if k in out: print('  %-16s %s' % (k, out[k]))
    return out
