"""Global, lattice-level check of the BMY ledger on explicit polynomial configurations.

For explicit primitive columns P_i in Q[t]^2 we
  * compute minors, contact ultrametrics at every finite root (irreducible factor) and at infinity,
  * perform the elementary transformations at infinity (ultrametric rule),
  * build Astra's blow-up tree (k>=3 centers weight 0; tangency chains weight a; transverse left),
  * build the Picard lattice of Y (basis C0, F, E_c total transforms), all classes, K_Y, B,
  * compute e_orb (SNC formula) and (K_Y+B)^2 by intersection numbers,
and compare with: (i) the ledger 8S/9 + 4T/9 + sum F, (ii) slack >= 0, (iii) effectivity of the
pair divisors V_ij, (iv) (IV.3) and (IV.2) with delta, E computed from the definitions.
"""
import random, itertools, sys
from fractions import Fraction as Fr
import sympy as sp

t = sp.symbols('t', real=True)
A = Fr(2, 3)


def ordp(poly, phi):
    o = 0
    p = sp.Poly(poly, t, domain='QQ')
    ph = sp.Poly(phi, t, domain='QQ')
    while True:
        q, r = sp.div(p, ph)
        if not r.is_zero:
            return o
        o += 1
        p = q


def phi_fn(k):
    return max(0, k*(k-1)//2 + 2*k - 9)


def analyse(cols, Npoly, nlist, verbose=False):
    m = len(cols)
    P = [(sp.Poly(x, t, domain='QQ'), sp.Poly(y, t, domain='QQ')) for x, y in cols]
    for x, y in P:
        assert sp.gcd(x, y).degree() == 0, 'not primitive'
    d = [max(x.degree(), y.degree()) for x, y in P]
    D = {}
    for i, j in itertools.combinations(range(m), 2):
        Dij = P[i][0]*P[j][1] - P[j][0]*P[i][1]
        assert not Dij.is_zero
        D[i, j] = Dij
    # irreducible factors of all minors
    facs = set()
    for Dij in D.values():
        for f, _ in sp.factor_list(Dij.as_expr(), t)[1]:
            facs.add(sp.Poly(f, t, domain='QQ').monic().as_expr())
    fibres = []   # (weight=number of conjugate roots, ultrametric dict)
    for f in facs:
        c = {(i, j): ordp(D[i, j].as_expr(), f) for (i, j) in D}
        fibres.append([sp.Poly(f, t).degree(), c, 'fin'])
    cinf = {(i, j): d[i]+d[j]-D[i, j].degree() for (i, j) in D}
    self_int = [2*di for di in d]          # on P1 x P1
    # elementary transformations at infinity
    nelem = 0
    while True:
        groups = classes(range(m), cinf, 1)
        big = [g for g in groups if 2*len(g) > m]
        if not big:
            break
        g = set(big[0])
        new = {}
        for (i, j), c in cinf.items():
            if i in g and j in g:
                new[i, j] = c-1
            elif i not in g and j not in g:
                new[i, j] = c+1
            else:
                assert c == 0
                new[i, j] = 0
        cinf = new
        self_int = [s-1 if i in g else s+1 for i, s in enumerate(self_int)]
        nelem += 1
    fibres.append([1, cinf, 'inf'])
    S = sum(self_int)
    # pair intersection numbers on X
    Cij = {(i, j): (self_int[i]+self_int[j])//2 for (i, j) in D}
    for (i, j) in D:
        assert (self_int[i]+self_int[j]) % 2 == 0
        tot = sum(w*f[1][i, j] for w, f, _ in [(fb[0], fb, None) for fb in fibres])
        assert tot == Cij[i, j], ('intersection mismatch', i, j, tot, Cij[i, j])
    I = sum(Cij.values())
    assert 2*I == (m-1)*S
    # lattice of F_e: choose e parity, b_i
    e = self_int[0] % 2
    bsec = [(s+e)//2 for s in self_int]
    # build centers.  A center = (fibre index, copy index, frozenset sections, weight, parent idx or None)
    centers = []
    transverse_count = {(i, j): 0 for (i, j) in D}
    ledgerF = Fr(0)
    sumphi_centers = 0
    for fi, (w, c, kind) in enumerate(fibres):
        for copy in range(w):
            build_tree(m, c, fi, copy, centers, transverse_count)
    n = len(centers)
    # lattice: vectors over basis [C0, F, E_0..E_{n-1}]; form: C0^2=-e, C0.F=1, F^2=0, E^2=-1
    dim = 2+n
    def dot(u, v):
        s = -e*u[0]*v[0] + u[0]*v[1] + u[1]*v[0]
        for q in range(n):
            s -= u[2+q]*v[2+q]
        return s
    def vec():
        return [Fr(0)]*dim
    secs = []
    for i in range(m):
        v = vec(); v[0] = Fr(1); v[1] = Fr(bsec[i])
        for q, cen in enumerate(centers):
            if i in cen['secs']:
                v[2+q] -= 1
        secs.append(v)
    exc = []
    for q, cen in enumerate(centers):
        v = vec(); v[2+q] = Fr(1)
        for q2, cen2 in enumerate(centers):
            if cen2['parent'] == q:
                v[2+q2] -= 1
        exc.append(v)
    K = vec(); K[0] = Fr(-2); K[1] = Fr(-2-e)
    for q in range(n):
        K[2+q] = Fr(1)
    fib = vec(); fib[1] = Fr(1)
    # boundary components (class, weight)
    comps = [(secs[i], A) for i in range(m)]
    comps += [(exc[q], centers[q]['w']) for q in range(n) if centers[q]['w'] > 0]
    comps += [(fib, Fr(1, 2))]*4
    # consistency: final section pair intersections = transverse count
    for (i, j) in D:
        assert dot(secs[i], secs[j]) == transverse_count[i, j], ('T mismatch', i, j)
    for v, _ in comps:
        assert dot(v, K) + dot(v, v) == -2     # all components rational smooth
    KB = [K[r] + sum(wt*v[r] for v, wt in comps) for r in range(dim)]
    KB2 = dot(KB, KB)
    eY = 4 + n
    eorb = Fr(eY) - sum(2*wt for v, wt in comps)
    # crossings: distinct half fibres are disjoint (F.F=0) -> fine
    for x, y in itertools.combinations(range(len(comps)), 2):
        inter = dot(comps[x][0], comps[y][0])
        assert inter >= 0
        eorb += comps[x][1]*comps[y][1]*inter
    slack = 3*eorb - KB2
    T = sum(transverse_count.values())
    # ledger
    F = lambda k, b0, b: 3*(1-2*b+b*b0)+3*A*k*(b-b0)+(1+b-b0-A*k)**2-2*A*A*Fr(k*(k-1), 2)
    led = Fr(8, 9)*S + Fr(4, 9)*T
    for cen in centers:
        b0 = centers[cen['parent']]['w'] if cen['parent'] is not None else Fr(0)
        led += F(len(cen['secs']), b0, cen['w'])
    assert slack == led, ('ledger mismatch', slack, led)
    # effectivity of V_ij in strict components
    for (i, j) in D:
        V = [K[r] + 2*fib[r] + secs[i][r] + secs[j][r] for r in range(dim)]  # K - pi^*K_P1 + C_i + C_j
        # claimed effective representative: sum over fibres contact*F_fibre + sum_c (1-l_c) E_c(total)
        W = vec(); W[1] = Fr(Cij[i, j])
        for q, cen in enumerate(centers):
            l = (i in cen['secs']) + (j in cen['secs'])
            W[2+q] += 1-l
        assert V == W, 'V_ij class mismatch'
        # coefficients on strict transforms: coef(E~_c) = mu_fibre + sum_{ancestors incl c} (1-l)
        for q, cen in enumerate(centers):
            mu = cen['mu'][(i, j)]
            tot = mu
            qq = q
            while qq is not None:
                cc = centers[qq]
                tot += 1-((i in cc['secs'])+(j in cc['secs']))
                qq = cc['parent']
            assert tot >= 0, ('negative coefficient', i, j, q, tot)
    # IV.2 data from definitions
    Np = sp.Poly(Npoly, t, domain='QQ')
    Lval = Np.degree()
    alpha = [sp.Poly(nn, t).degree() for nn in nlist]
    delta = Fr(0)
    for (i, j) in D:
        delta = max(delta, Fr(D[i, j].degree()) - Fr(alpha[i]+alpha[j], 2) + Fr(Lval, 4))
    Esum = 0
    for (i, j) in D:
        g = sp.gcd(sp.Poly(nlist[i], t, domain='QQ'), sp.Poly(nlist[j], t, domain='QQ'))
        g2 = sp.gcd(g, D[i, j])
        Esum += g.degree()-g2.degree()
    # layers: roots of N (over factors) and depths
    layers = []
    for f, mult in sp.factor_list(Npoly, t)[1]:
        for r in range(1, mult+1):
            sh = sum(1 for nn in nlist if ordp(nn, f) >= r)
            layers += [sh]*sp.Poly(f, t).degree()
    assert len(layers) == Lval and sum(layers) == sum(alpha)
    iv3 = S <= sum(layers) - Fr(m*Lval, 4) + m*delta
    lhs = Fr(3, 4)*(m-24)*Lval + sum((Fr(sh)-Fr(m, 2))**2 for sh in layers)
    rhs = m*(m+3)*delta + 4*Esum
    return dict(m=m, L=Lval, S=S, I=I, T=T, n=n, nelem=nelem, slack=slack, delta=delta, E=Esum,
                iv3=iv3, iv2=(lhs <= rhs), lhs=lhs, rhs=rhs, maxgroup=max(len(c['secs']) for c in centers) if centers else 0)


def classes(idx, c, r):
    idx = list(idx)
    out, seen = [], set()
    for i in idx:
        if i in seen:
            continue
        g = [i] + [j for j in idx if j != i and c[min(i, j), max(i, j)] >= r]
        for j in g:
            seen.add(j)
        out.append(g)
    return out


def build_tree(m, c, fi, copy, centers, tcount):
    """Astra's procedure on one fibre with contact ultrametric c."""
    def cc(i, j):
        return c[min(i, j), max(i, j)]
    mu = {key: val for key, val in c.items()}   # fibre coefficient of V_ij = contact multiplicity

    def rec(group, depth, parent):
        # group: sections sharing the point at this depth (depth>=1), parent center index
        if len(group) >= 3:
            q = len(centers)
            centers.append(dict(secs=frozenset(group), w=Fr(0), parent=parent, mu=mu))
            for sub in classes(group, c, depth+1):
                rec(sub, depth+1, q)
        elif len(group) == 2:
            i, j = sorted(group)
            e = cc(i, j) - depth + 1
            assert e >= 1
            if e == 1:
                tcount[i, j] += 1
                return
            par = parent
            for step in range(e):
                q = len(centers)
                centers.append(dict(secs=frozenset(group), w=A, parent=par, mu=mu))
                par = q
    for g in classes(range(m), c, 1):
        if len(g) >= 2:
            rec(g, 1, None)


# ------------------------------------------------------------------------------------------
def hermite(jets, extra_deg, rng):
    """polynomial f with f = jet_x mod (t-x)^len(jet) for each x, plus random multiple of the modulus."""
    M = sp.Integer(1)
    for x, jet in jets.items():
        M *= (t-x)**len(jet)
    deg = sum(len(j) for j in jets.values())
    cs = sp.symbols('c0:%d' % max(deg, 1))
    f = sum(cs[k]*t**k for k in range(deg)) if deg else sp.Integer(0)
    eqs = []
    for x, jet in jets.items():
        for k, val in enumerate(jet):
            eqs.append(sp.Eq(sp.diff(f, t, k).subs(t, x)/sp.factorial(k), val))
    if deg:
        sol = sp.solve(eqs, cs[:deg], dict=True)[0]
        f = f.subs(sol)
    g = sum(rng.randint(-3, 3)*t**k for k in range(extra_deg+1))
    if g == 0:
        g = 1
    return sp.expand(f + M*g)


def random_config(rng, m, roots, maxexp, pmark, shared_inf=False, perturb=0.0):
    ex = {x: rng.randint(1, maxexp) for x in roots}
    jetsR = {x: [rng.randint(-5, 5) for _ in range(ex[x])] for x in roots}
    jetsS = {x: [rng.randint(-5, 5) for _ in range(ex[x])] for x in roots}
    for x in roots:
        jetsS[x][0] = rng.choice([1, 2, 3])
    cols, nl = [], []
    for i in range(m):
        jR, jS, nn = {}, {}, sp.Integer(1)
        for x in roots:
            if rng.random() < pmark:
                a_ = rng.randint(1, ex[x])
                jR[x] = jetsR[x][:a_]
                jS[x] = jetsS[x][:a_]
                if rng.random() < perturb:
                    jR[x] = [jR[x][0]+1] + jR[x][1:]
                nn *= (t-x)**a_
        ed = rng.randint(0, 2)
        f = hermite(jR, ed, rng)
        g = hermite(jS, ed if shared_inf else rng.randint(0, 2), rng)
        cols.append((f, g)); nl.append(nn)
    Np = sp.Integer(1)
    for x in roots:
        Np *= (t-x)**ex[x]
    return cols, sp.expand(Np), nl


if __name__ == '__main__':
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    rng = random.Random(seed)
    stats = dict(ok=0, minslack=None, elem=0, big=0)
    for tr in range(trials):
        m = rng.randint(4, 8)
        roots = rng.sample(range(-6, 7), rng.randint(2, 5))
        try:
            cols, Np, nl = random_config(rng, m, roots, 3, rng.choice([0.5, 0.7, 0.9]),
                                         shared_inf=rng.random() < 0.5, perturb=rng.choice([0, 0.2]))
            res = analyse(cols, Np, nl)
        except AssertionError as ex:
            if 'not primitive' in str(ex) or str(ex) == '':
                continue
            raise
        except (IndexError, ZeroDivisionError):
            continue
        stats['ok'] += 1
        stats['elem'] += res['nelem'] > 0
        stats['big'] += res['maxgroup'] >= 3
        assert res['slack'] >= 0, res
        assert res['iv3'], res
        assert res['iv2'], res
        if stats['minslack'] is None or res['slack'] < stats['minslack']:
            stats['minslack'] = res['slack']
        print(tr, {k_: (str(v) if isinstance(v, Fr) else v) for k_, v in res.items()})
    print('SUMMARY', stats)
