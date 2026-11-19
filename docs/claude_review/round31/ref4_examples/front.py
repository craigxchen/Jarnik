"""sympy front end: exact contact data and (IV.2) data from polynomial columns."""
import sympy as sp
from engine import chern_side, iv2, classes
from fractions import Fraction as Fr

t = sp.symbols('t')


def P(e, ext):
    if ext is None:
        return sp.Poly(e, t, domain='QQ')
    return sp.Poly(e, t, extension=ext)


def vval(D, p):
    v = 0
    q = D
    while True:
        qq, r = sp.div(q, p)
        if not r.is_zero:
            return v
        v += 1
        q = qq


def key(p):
    return str(p.monic().as_expr())


def contact_data(cols, ext=None):
    m = len(cols)
    X = [(P(x, ext), P(y, ext)) for x, y in cols]
    for x, y in X:
        assert sp.gcd(x, y).degree() == 0, "column not primitive"
    degs = [max(x.degree(), y.degree()) for x, y in X]
    D = {}
    facs = {}
    for i in range(m):
        for j in range(i + 1, m):
            d = X[i][0] * X[j][1] - X[j][0] * X[i][1]
            assert not d.is_zero, ("zero minor", i, j)
            D[(i, j)] = d
            for f, e in d.factor_list()[1]:
                if f.degree() > 0:
                    facs.setdefault(key(f), f.monic())
    finite = []
    vals = {}
    for k, p in facs.items():
        v = [[0] * m for _ in range(m)]
        for (i, j), d in D.items():
            v[i][j] = v[j][i] = vval(d, p)
        vals[k] = v
        finite.append((p.degree(), v))
    vinf = [[0] * m for _ in range(m)]
    for (i, j), d in D.items():
        vinf[i][j] = vinf[j][i] = degs[i] + degs[j] - d.degree()
    return dict(m=m, X=X, D=D, degs=degs, facs=facs, vals=vals, finite=finite, vinf=vinf)


def greedy_marks(cd, keys=None, minsize=2):
    """For each finite contact factor p choose the nested chain of largest
    classes; returns N = {key: (p, e_p)}, a = list of {key: a_ip}."""
    m = cd['m']
    N = {}
    a = [dict() for _ in range(m)]
    for k, p in cd['facs'].items():
        if keys is not None and k not in keys:
            continue
        v = cd['vals'][k]
        group = list(range(m))
        d = 0
        while True:
            cls = [c for c in classes(group, v, d + 1) if len(c) >= minsize]
            if not cls:
                break
            g = max(cls, key=len)
            d += 1
            for i in g:
                a[i][k] = d
            group = g
        if d > 0:
            N[k] = (p, d)
    return N, a


def iv2_data(cd, N, a, delta=None):
    """N: {key: (p, e)}; a: list of {key: exponent}.  p must be irreducible
    over the base field (so that its roots form one layer orbit)."""
    m = cd['m']
    L = sum(p.degree() * e for p, e in N.values())
    layers = []
    for k, (p, e) in N.items():
        for d in range(1, e + 1):
            s = sum(1 for i in range(m) if a[i].get(k, 0) >= d)
            layers.append((p.degree(), s))
    alpha = [sum(N[k][0].degree() * x for k, x in a[i].items()) for i in range(m)]
    for i in range(m):
        for k, x in a[i].items():
            assert x <= N[k][1]
    degD = {}
    E = {}
    for (i, j), d in cd['D'].items():
        degD[(i, j)] = d.degree()
        tot = 0
        for k, (p, e) in N.items():
            g = min(a[i].get(k, 0), a[j].get(k, 0))
            if g:
                v = cd['vals'][k][i][j] if k in cd['vals'] else vval(d, p)
                tot += p.degree() * (g - min(g, v))
        E[(i, j)] = tot
    return iv2(m, L, layers, alpha, degD, E, delta)


def full(cols, ext=None, N=None, a=None, label=''):
    cd = contact_data(cols, ext)
    ch = chern_side(cd['m'], cd['finite'], cd['vinf'], cd['degs'])
    if N is None:
        N, a = greedy_marks(cd)
    ar = iv2_data(cd, N, a)
    return cd, ch, ar


def report(label, ch, ar):
    s = (f"{label}: m={len(ar) and ''}"
         f"S={ch['S']} I={ch['I']} infsteps={ch['steps']} blowups={ch['nblow']} "
         f"big={ch['big'][:8]}{'...' if len(ch['big'])>8 else ''} chains={ch['chains'][:8]} unblown={ch['unblown']} | "
         f"(K+B)^2={ch['KB2']} e_orb={ch['eorb']} 3e-(K+B)^2={ch['Q']} ledger={ch['ledger']} bound={ch['bound']} "
         f"min(K+B).curve={ch['nefmin']} || IV.2: L={ar['L']} delta={ar['delta']} sumE={ar['sumE']} "
         f"lhs={ar['lhs']} rhs={ar['rhs']} slack={ar['slack']}")
    print(s, flush=True)
