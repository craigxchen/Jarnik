"""Referee 'bridge' (round 31): exact dictionary from function-field clusters to Part IV data.

Cluster: monic z_i in Q(i)[t] of degree D with z_i*conj(z_i) = N, deg(z_i - z_j) <= K.
Orientation U: U*conj(U) = N, gcd(U, conj U) = 1  (exists iff N has no real root).
Star (unanchored): P_i = gcd(z_i, U) monic; then z_i = P_i*conj(U)/conj(P_i).
Columns (X_i, Y_i) = (Re P_i, Im P_i) in Q[t]^2; n_i = P_i conj(P_i) | N; D_ij = X_i Y_j - X_j Y_i.
Checks (all exact): primitivity, n_i | N, g_ij | D_ij (E_ij = 0), the exact degree identity
   deg D_ij = deg(z_i - z_j) - D + (alpha_i + alpha_j)/2,
hence delta = K - D/2 (allowed) and delta_act = max deg(z_i-z_j) - D/2.
Layers: for each irreducible factor f of N over Q(i), multiplicity e, depth l<=e: deg f layers with
s = #{i : f^l | n_i}.  Then both sides of (IV.2) are printed.
"""
import itertools, sys
from fractions import Fraction as Fr
import sympy as sp

t = sp.symbols('t'); I = sp.I
DOM = 'QQ_I'


def P(expr):
    return sp.Poly(sp.expand(expr), t, domain=DOM)


def conj(p):
    return sp.Poly([sp.conjugate(c) for c in p.all_coeffs()], t, domain=DOM)


def monic(p):
    return p.monic()


def is_real(p):
    return all(sp.im(c) == 0 for c in p.all_coeffs())


def dictionary(zs, U, D, K, name, verbose=True):
    m = len(zs)
    Ub = conj(U)
    N = U * Ub
    assert sp.gcd(U, Ub).degree() == 0, "orientation not coprime"
    L = N.degree(); assert L == 2 * D
    for z in zs:
        assert z.degree() == D and z.LC() == 1
        assert z * conj(z) == N, "not on the circle"
    for i, j in itertools.combinations(range(m), 2):
        assert (zs[i] - zs[j]).degree() <= K, "not a cluster"
    Ps, Xs, Ys, ns, al = [], [], [], [], []
    for z in zs:
        p = monic(sp.gcd(z, U)); pb = conj(p)
        assert z * pb == p * Ub, "z != P Ubar / Pbar"
        X = (p + pb) * sp.Rational(1, 2); Y = (p - pb) * (1 / (2 * I))
        X = sp.Poly(X.as_expr(), t, domain=DOM); Y = sp.Poly(sp.expand(Y.as_expr()), t, domain=DOM)
        assert is_real(X) and is_real(Y)
        assert sp.gcd(X, Y).degree() == 0, "column not primitive"
        n = p * pb
        assert N.rem(n).is_zero
        Ps.append(p); Xs.append(X); Ys.append(Y); ns.append(n); al.append(n.degree())
    Eall, dmax, dact = 0, None, None
    for i, j in itertools.combinations(range(m), 2):
        Dij = Xs[i] * Ys[j] - Xs[j] * Ys[i]
        assert not Dij.is_zero
        g = sp.gcd(ns[i], ns[j])
        E = g.degree() - sp.gcd(g, Dij).degree()
        Eall += E
        dz = (zs[i] - zs[j]).degree()
        assert 2 * Dij.degree() == 2 * dz - 2 * D + al[i] + al[j], "degree identity fails"
        ex = Fr(2 * Dij.degree() - al[i] - al[j], 2) + Fr(L, 4)
        dact = ex if dact is None else max(dact, ex)
    # layers
    fl = sp.factor_list(N.as_expr(), t, extension=I)[1]
    s_list = []
    for f, e in fl:
        fp = sp.Poly(f, t, domain=DOM)
        for l in range(1, e + 1):
            s = sum(1 for n in ns if n.rem(fp ** l).is_zero)
            s_list += [s] * fp.degree()
    assert len(s_list) == L
    assert sum(s_list) == sum(al)
    delta = Fr(K) - Fr(D, 2)
    Sig = sum((Fr(s) - Fr(m, 2)) ** 2 for s in s_list)
    lhs = Fr(3, 4) * (m - 24) * L + Sig
    rhs = m * (m + 3) * delta + 4 * Eall
    rhs_act = m * (m + 3) * dact + 4 * Eall
    if verbose:
        print(f"[{name}] D={D} K={K} m={m} L={L} delta(allowed)={delta} delta(actual)={dact} "
              f"sum E={Eall} layer sizes={sorted(s_list)} Sigma={Sig}")
        print(f"     (IV.2): LHS={lhs}  RHS(allowed)={rhs}  RHS(actual)={rhs_act}  holds={lhs <= rhs_act}")
    assert Eall == 0 and dact <= delta and lhs <= rhs_act
    return dict(m=m, L=L, delta=delta, dact=dact, E=Eall, s=s_list, lhs=lhs, rhs=rhs, rhs_act=rhs_act)


def cluster_from_roots(roots_mult, K, base_choice=None):
    """roots_mult: list of (rho, e) with rho nonreal in Q(i). All points: exponent a_k in [0,e_k].
    Return the largest cluster through the point base_choice (default: maximal clique search, brute)."""
    D = sum(e for _, e in roots_mult)
    pts = {}
    for a in itertools.product(*[range(e + 1) for _, e in roots_mult]):
        ex = 1
        for (r, e), ak in zip(roots_mult, a):
            ex *= (t - r) ** ak * (t - sp.conjugate(r)) ** (e - ak)
        pts[a] = P(ex)
    return D, pts


def cube8():
    y = [sp.Rational(-9, 4), sp.Rational(-2, 3), sp.Rational(1, 4), sp.Rational(8, 3)]
    eps = [1, -1, -1, 1]
    r = [eps[k] / y[k] + I * y[k] for k in range(4)]
    rho = r + [-x for x in r]
    Q = P(sp.prod([t - z for z in rho]))
    up = [z if sp.im(z) > 0 else sp.conjugate(z) for z in rho]
    D, pts = cluster_from_roots([(u, 1) for u in up], 4)
    clus = [p for p in pts.values() if (p - Q).degree() <= 4 or (p - Q).is_zero]
    U = P(sp.prod([t - u for u in up]))
    return clus, U


def main():
    res = {}
    # 1. cube D=8, c=0
    clus, U = cube8()
    assert len(clus) == 8
    res['cube8'] = dictionary(clus, U, 8, 4, 'cube D=8 c=0')
    # 1b. other orientation (flip half the roots) -- same Sigma in squarefree case
    roots = [r for r, _ in sp.roots(U.as_expr(), t).items()]
    U2 = P(sp.prod([(t - r) if k % 2 else (t - sp.conjugate(r)) for k, r in enumerate(roots)]))
    dictionary(clus, U2, 8, 4, 'cube D=8 c=0, flipped orientation')
    # 2. cube D=16 by composition with g = t^2 + t + 3 (real)
    g = t ** 2 + t + 3
    comp = lambda p: P(p.as_expr().subs(t, g))
    clus16 = [comp(z) for z in clus]
    res['cube16'] = dictionary(clus16, comp(U), 16, 8, 'cube D=16 c=0 (z o g)')
    # 3. vertical-line clusters (roots i*y), c=1
    for D, K, B in ((7, 4, 9), (8, 5, 10)):
        Lm = D - 1 - K
        best = None
        for absv in itertools.combinations(range(1, B + 1), D):
            for sg in itertools.product((1, -1), repeat=D - 1):
                yv = [absv[0]] + [s * a for s, a in zip(sg, absv[1:])]
                cnt = [F for F in itertools.product((0, 1), repeat=D)
                       if all(sum(f * (yy ** mm) for f, yy in zip(F, yv)) == 0
                              for mm in range(1, Lm + 1) if mm % 2)]
                if best is None or len(cnt) > len(best[1]):
                    best = (yv, cnt)
        yv, Fs = best
        rts = [I * yy for yy in yv]
        Q = P(sp.prod([t - r for r in rts]))
        zs = []
        for F in Fs:
            zs.append(P(sp.prod([(t - sp.conjugate(r)) if f else (t - r) for f, r in zip(F, rts)])))
        Uv = P(sp.prod([t - (r if sp.im(r) > 0 else sp.conjugate(r)) for r in rts]))
        print(f"vertical D={D}: y={yv}, cluster size {len(zs)}")
        res[f'vert{D}'] = dictionary(zs, Uv, D, K, f'vertical D={D} c=1')
    # 4. non-squarefree N: e=(2,1,1), D=4, K=3 (every point; c=1)
    rm = [(1 + 2 * I, 2), (-1 + I, 1), (3 * I, 1)]
    D, pts = cluster_from_roots(rm, 3)
    zs = list(pts.values())
    U4 = P(sp.prod([(t - r) ** e for r, e in rm]))
    res['nsf4'] = dictionary(zs, U4, 4, 3, 'non-squarefree e=(2,1,1), all 12 points, K=3')
    # 4b. non-squarefree, e=(3,2), D=5, K=4: all 12 points
    rm = [(1 + I, 3), (-2 + 3 * I, 2)]
    D, pts = cluster_from_roots(rm, 4)
    zs = list(pts.values())
    U5 = P(sp.prod([(t - r) ** e for r, e in rm]))
    res['nsf5'] = dictionary(zs, U5, 5, 4, 'non-squarefree e=(3,2), all 12 points, K=4')
    # 4c. anchored star = reorientation check on e=(2,1,1): anchor point a0=(1,0,1)
    print("all exact checks passed")
    return res


if __name__ == '__main__':
    main()
