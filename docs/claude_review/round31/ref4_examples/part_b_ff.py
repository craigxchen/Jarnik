"""(b) The record's function-field clusters passed to star matrices.
Model (round29/30): N real squarefree of degree 2D, z monic of degree D in C[t]
with z zbar = N, cluster condition deg(z_i - z_j) <= K = floor(D/2) + c.
Star of z_F: P_F = prod of the base roots kept; column (Re P_F, Im P_F)."""
import itertools, sys
import sympy as sp
from front import t, contact_data, iv2_data, key, report
from engine import chern_side

I = sp.I


def star_config(base_roots, sets, D, K, h=None):
    """base_roots: list of D roots rho_k (all non-real, pairwise distinct with conjugates);
    sets: list of flip sets F (subsets of range(D)).  Returns columns, N, a."""
    z = {}
    cols = []
    for F in sets:
        keep = [k for k in range(D) if k not in F]
        Pp = sp.expand(sp.prod([t - base_roots[k] for k in keep]))
        zz = sp.expand(sp.prod([t - (sp.conjugate(base_roots[k]) if k in F else base_roots[k]) for k in range(D)]))
        if h is not None:
            Pp = sp.expand(Pp.subs(t, h))
            zz = sp.expand(zz.subs(t, h))
        z[F] = zz
        re = sp.expand((Pp + sp.conjugate(Pp).subs(sp.conjugate(t), t)) / 2)
        im = sp.expand((Pp - sp.conjugate(Pp).subs(sp.conjugate(t), t)) / (2 * I))
        cols.append((re, im))
    # verify the cluster condition exactly
    mx = 0
    for F, G in itertools.combinations(sets, 2):
        dg = sp.Poly(sp.expand(z[F] - z[G]), t).degree()
        mx = max(mx, dg)
    Dd = D * (1 if h is None else sp.Poly(h, t).degree())
    assert mx <= K, (mx, K)
    # primes p_k = (t-rho)(t-rhobar) (composed with h), factored over Q
    N = {}
    a = [dict() for _ in sets]
    for k in range(D):
        pk = sp.expand((t - base_roots[k]) * (t - sp.conjugate(base_roots[k])))
        if h is not None:
            pk = sp.expand(pk.subs(t, h))
        for f, e in sp.Poly(pk, t, domain='QQ').factor_list()[1]:
            kk = key(f)
            N[kk] = (f.monic(), e)
            for i, F in enumerate(sets):
                if k not in F:
                    a[i][kk] = e
    return cols, N, a, mx, Dd


def go(label, base_roots, sets, D, K, h=None):
    cols, N, a, mx, Dd = star_config(base_roots, sets, D, K, h)
    cd = contact_data(cols)
    ar = iv2_data(cd, N, a)
    ch = chern_side(cd['m'], cd['finite'], cd['vinf'], cd['degs'])
    print(f"--- {label}: D={Dd} K={K} max deg(z_i-z_j)={mx} m={len(sets)} predicted delta=K-D/2={sp.Rational(K) - sp.Rational(Dd, 2)}")
    report(label, ch, ar)
    return ch, ar


if __name__ == '__main__':
    # cube cluster D = 8, c = 0 (round 29 Theorem C), exact over Q(i)
    y = [sp.Rational(-9, 4), sp.Rational(-2, 3), sp.Rational(1, 4), sp.Rational(8, 3)]
    eps = [1, -1, -1, 1]
    r = [eps[k] / y[k] + I * y[k] for k in range(4)]
    roots = r + [-x for x in r]          # roots of Q; index k and k+4 = -r_k
    pm = lambda ks: frozenset(ks + [k + 4 for k in ks])
    sets = [frozenset(), frozenset(range(8)), pm([0, 1]), pm([2, 3]), pm([0, 2]), pm([1, 3]),
            frozenset([0, 1, 2, 3]), frozenset([4, 5, 6, 7])]
    go("cube D=8", roots, sets, 8, 4)
    # D = 16 by composition with h = t^2 + t
    go("cube D=16 (h=t^2+t)", roots, sets, 8, 8, h=t ** 2 + t)

    # c = 1 vertical-line clusters: rho_k = i*Y_k; flip set F zero-sum iff sum_F Y_k = 0 (L = 2)
    D = 8
    best = None
    for absv in itertools.combinations(range(1, 11), D):
        for sg in itertools.product([1, -1], repeat=D - 1):
            Y = [absv[0]] + [s * x for s, x in zip(sg, absv[1:])]
            cnt = sum(1 for F in itertools.product([0, 1], repeat=D) if sum(f * yy for f, yy in zip(F, Y)) == 0)
            if best is None or cnt > best[0]:
                best = (cnt, Y)
    print("best vertical-line Y at D=8, |y|<=10:", best)
    Y = best[1]
    sets = [frozenset(k for k in range(D) if F[k]) for F in itertools.product([0, 1], repeat=D)
            if sum(f * yy for f, yy in zip(F, Y)) == 0]
    go("vertical c=1 D=8", [I * yy for yy in Y], sets, 8, 5)
