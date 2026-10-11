"""Exact inner analysis of a weighted support (T, om) for the normalised problem.

analyse(M, T, om, rho0, rho1, strict) enumerates, in exact arithmetic, every class of G_T (integral
class included) and every lift lambda with sigma_lambda < r(lambda) (or <= if strict), using the
coercivity box of verify.py.  It returns
  * 'rank_violator' if the rows are affinely dependent on T;
  * 'int_viol': violated integral characters;
  * 'bad': {element n of G_T: (class key, [violated lifts])} for non-integral classes.
Then find_hom(G, bad) searches a homomorphism h: G_T -> Q/Z with h(n) != 0 for every bad n (windings
realising it: G.kappa_of(eta)); none exists iff the bad classes form a *covering family*.

Lemma (covering families are intrinsic): for a finite family F of nonzero classes in Q^M_0/Z^M_0 and
any support T' with F inside G_T', "every homomorphism G_T' -> Q/Z vanishes on some member of F" holds
iff it holds for the subgroup <F> (homomorphisms of <F> extend to G_T' because Q/Z is injective).
So a covering family found for one support is covering for every support whose lattice contains it.
"""
import math, itertools
from fractions import Fraction
import numpy as np
from core import Group, Qmat, frac_part, left_inverse_lp, exact_left_inverse, coercivity, enum_class, lcm_list
from verify import integer_kernel_vector


def analyse(M, T, om, rho0=0, rho1=0, strict=False, max_lifts=6, max_int=300, max_group=200000):
    T = list(T)
    om = [Fraction(w) for w in om]
    rho0 = Fraction(rho0); rho1 = Fraction(rho1)
    Om = sum(om)
    Q = Qmat(M, T)
    out = {'Omega': Om}
    G = Group(M, T)
    out['G'] = G
    if not G.full_rank:
        out['rank_violator'] = integer_kernel_vector(Q, M - 1)
        return out
    if G.order > max_group:
        out['too_large'] = True
        return out
    A = exact_left_inverse(Q, left_inverse_lp(Q, [float(w) for w in om]))
    k_inf, k_1 = coercivity(Q, om, A)
    if rho1 > 0 and (rho1 * k_1 >= Fraction(3, 2) or k_1 > 1):
        from nu1 import kappa1_bound
        kb = kappa1_bound(M, T, om)
        if kb is not None and kb < k_1:
            k_1 = kb
    k_inf = min(k_inf, k_1 / 2)        # ||lambda||_inf <= ||lambda||_1 / 2 for zero-sum lambda
    if rho1 * k_1 >= 2:
        out['inconclusive'] = True
        return out
    Nmax = (Om / 2 + rho0) / (1 - rho1 * k_1 / 2)
    Binf = k_inf * Nmax
    B1 = k_1 * Nmax
    E = lcm_list([w.denominator for w in om])
    Wint = np.array([int(w * E) for w in om], dtype=np.int64)
    SW = int(Wint.sum())
    F0 = rho0.denominator; a0 = rho0.numerator
    F1 = rho1.denominator; a1 = rho1.numerator
    QT = Q.T.astype(np.int64)
    int_viol = []
    bad = {}
    min_slack = None
    elems = list(G.elements())
    for nv in elems:
        if any(nv):
            mu = G.element(nv)
            f = frac_part(mu)
            D = lcm_list([v.denominator for v in f])
            F = [int(v * D) for v in f]
        else:
            f = [Fraction(0)] * (M - 1); D = 1; F = [0] * (M - 1)
        scale = 2 * E * D * F0 * F1
        lifts = []
        for V in enum_class(F, D, math.floor(Binf * D), math.floor(B1 * D)):
            LS = np.abs(V @ QT)
            l1 = np.abs(V).sum(axis=1) + np.abs(V.sum(axis=1))
            val = (2 * F0 * F1) * (LS @ Wint) - F0 * F1 * D * SW - 2 * E * D * a0 * F1 - E * a1 * F0 * l1
            j = int(np.argmin(val))
            m = Fraction(int(val[j]), scale)
            if min_slack is None or m < min_slack[0]:
                min_slack = (m, [Fraction(int(v), D) for v in V[j]])
            mask = (val <= 0) if strict else (val < 0)
            if mask.any():
                idx = np.where(mask)[0]
                idx = idx[np.argsort(val[idx])][:max_lifts if any(nv) else max_int]
                for jj in idx:
                    lifts.append((Fraction(int(val[jj]), scale), [Fraction(int(v), D) for v in V[jj]]))
        if lifts:
            if any(nv):
                bad[tuple(nv)] = (tuple(f), lifts)
            else:
                int_viol.extend(lifts)
    out['int_viol'] = int_viol
    out['bad'] = bad
    out['min_slack'] = min_slack
    return out


def hom_matrix(G, elems):
    """values L*h_eta(n) mod L for all eta (rows) and n in elems (cols); L = lcm(invs)."""
    invs = G.invs
    L = lcm_list(invs) if invs else 1
    etas = np.array(list(itertools.product(*[range(d) for d in invs])), dtype=np.int64).reshape(-1, len(invs))
    N = np.array(elems, dtype=np.int64).reshape(-1, len(invs))
    coef = np.array([L // d for d in invs], dtype=np.int64)
    vals = (etas * coef) @ N.T % L
    return etas, vals


def find_hom(G, bad_elems):
    if not bad_elems:
        return [0] * len(G.invs)
    etas, vals = hom_matrix(G, bad_elems)
    good = np.all(vals != 0, axis=1)
    if good.any():
        return [int(v) for v in etas[np.argmax(good)]]
    return None


def minimal_cover(G, bad_elems):
    """Greedy inclusion-minimal covering subfamily (assumes bad_elems is covering)."""
    etas, vals = hom_matrix(G, bad_elems)
    keep = list(range(len(bad_elems)))
    zero = (vals == 0)
    for i in sorted(range(len(bad_elems)), key=lambda i: -zero[:, i].sum()):
        trial = [j for j in keep if j != i]
        if not trial:
            continue
        if np.all(zero[:, trial].any(axis=1)):
            keep = trial
    return [bad_elems[j] for j in keep]
