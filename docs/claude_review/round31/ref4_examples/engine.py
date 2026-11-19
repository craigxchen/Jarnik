"""Independent engine for refereeing Part IV, inequality (IV.2).

Two independent computations are done for a configuration of m sections of a
ruled surface over P^1 (columns P_i in K[t]^2, primitive):

 (1) the Chern-number side: Astra's surface Y (elementary transformations at
     infinity until no point of the infinity fibre carries > m/2 sections; blow up
     every centre with k>=3 sections, weight 0; blow up whole chains of pair
     tangencies of remaining order e>=2, weight a; four generic half fibres),
     with (K_Y+B)^2 and e_orb(Y,B) computed DIRECTLY from an intersection form on
     the span of pullback classes and total transforms of exceptional curves, and
     strict transforms written out explicitly.  This is compared with Astra's
     telescoped ledger and with the bound used in the proof.
 (2) both sides of (IV.2) from the arithmetic data (N, n_i, D_ij).

Exact rational arithmetic throughout (fractions.Fraction).
"""
from fractions import Fraction as Fr
from math import comb
import itertools

A = Fr(2, 3)


def phi(k):
    return max(0, comb(k, 2) + 2 * k - 9)


def check_ultrametric(v, m):
    for i, j, k in itertools.permutations(range(m), 3):
        assert v[i][k] >= min(v[i][j], v[j][k]), ("not ultrametric", i, j, k)


def classes(group, v, d):
    """partition `group` by the relation v[i][j] >= d (an equivalence relation
    for an ultrametric contact matrix)."""
    out = []
    left = list(group)
    while left:
        i = left.pop(0)
        cl = [i] + [j for j in left if v[i][j] >= d]
        left = [j for j in left if j not in cl]
        out.append(cl)
    return out


def normalise_infinity(m, vinf, selfint):
    """Elementary transformations at infinity while some infinity point carries
    k > m/2 sections.  Returns new vinf, selfint, number of transformations."""
    v = [row[:] for row in vinf]
    c2 = list(selfint)
    steps = 0
    while True:
        big = [cl for cl in classes(range(m), v, 1) if 2 * len(cl) > m]
        if not big:
            return v, c2, steps
        cl = set(big[0])
        steps += 1
        for i in range(m):
            c2[i] += -1 if i in cl else 1
        for i in range(m):
            for j in range(m):
                if i == j:
                    continue
                if i in cl and j in cl:
                    v[i][j] -= 1
                elif i not in cl and j not in cl:
                    v[i][j] += 1
                else:
                    assert v[i][j] == 0
                    v[i][j] = 0


class Surface:
    def __init__(self, m):
        self.m = m
        self.blow = []      # list of dicts: secs (frozenset), b, parent (index or None)
        self.unblown = 0    # number of unblown transverse pair points
        self.unblown_pairs = {}
        self.big = []       # sizes of k>=3 centres
        self.chains = []    # lengths e of pair chains
        self.I = 0

    def process_point(self, v, weight=1):
        m = self.m
        for i in range(m):
            for j in range(i + 1, m):
                self.I += weight * v[i][j]
        for _ in range(weight):
            for cl in classes(range(m), v, 1):
                if len(cl) >= 2:
                    self._centre(cl, 1, None, Fr(0), v)

    def _centre(self, group, d, parent, b0, v):
        k = len(group)
        if k >= 3:
            self.blow.append(dict(secs=frozenset(group), b=Fr(0), parent=parent))
            me = len(self.blow) - 1
            self.big.append(k)
            for cl in classes(group, v, d + 1):
                if len(cl) >= 2:
                    self._centre(cl, d + 1, me, Fr(0), v)
        else:
            i, j = group
            e = v[i][j] - d + 1
            assert e >= 1
            if e == 1:
                self.unblown += 1
                key = (min(i, j), max(i, j))
                self.unblown_pairs[key] = self.unblown_pairs.get(key, 0) + 1
            else:
                self.chains.append(e)
                par, pb = parent, b0
                for _ in range(e):
                    self.blow.append(dict(secs=frozenset(group), b=A, parent=par))
                    par = len(self.blow) - 1
                    pb = A

    # ---------- direct intersection theory ----------
    def direct(self, c2, Cdot):
        """c2[i] = C_i^2, Cdot[i][j] = C_i.C_j on the (normalised) ruled surface."""
        m = self.m
        nb = len(self.blow)
        # basis: 0:K0, 1:F, 2..2+m-1: C_i, then E_k
        n = 2 + m + nb
        G = [[Fr(0)] * n for _ in range(n)]
        G[0][0] = Fr(8)
        G[0][1] = G[1][0] = Fr(-2)
        G[1][1] = Fr(0)
        for i in range(m):
            G[0][2 + i] = G[2 + i][0] = Fr(-2 - c2[i])
            G[1][2 + i] = G[2 + i][1] = Fr(1)
            for j in range(m):
                G[2 + i][2 + j] = Fr(c2[i]) if i == j else Fr(Cdot[i][j])
        for k in range(nb):
            G[2 + m + k][2 + m + k] = Fr(-1)

        def dot(u, w):
            s = Fr(0)
            for p, x in u.items():
                if x == 0:
                    continue
                row = G[p]
                for q, y in w.items():
                    if y:
                        s += x * y * row[q]
            return s

        def add(u, w, c=Fr(1)):
            r = dict(u)
            for q, y in w.items():
                r[q] = r.get(q, Fr(0)) + c * y
            return r

        Cs = []
        for i in range(m):
            d = {2 + i: Fr(1)}
            for k, bl in enumerate(self.blow):
                if i in bl['secs']:
                    d[2 + m + k] = d.get(2 + m + k, Fr(0)) - 1
            Cs.append(d)
        Es = []
        for k, bl in enumerate(self.blow):
            d = {2 + m + k: Fr(1)}
            for l, bl2 in enumerate(self.blow):
                if bl2['parent'] == k:
                    d[2 + m + l] = d.get(2 + m + l, Fr(0)) - 1
            Es.append(d)
        Fv = {1: Fr(1)}
        KY = {0: Fr(1)}
        for k in range(nb):
            KY[2 + m + k] = Fr(1)
        comps = [(Cs[i], A) for i in range(m)]
        comps += [(Es[k], bl['b']) for k, bl in enumerate(self.blow) if bl['b'] > 0]
        comps += [(Fv, Fr(1, 2))] * 4
        Bd = {}
        for D, b in comps:
            Bd = add(Bd, D, b)
        KB = add(KY, Bd)
        KB2 = dot(KB, KB)
        e = Fr(4 + nb)
        for D, b in comps:
            e -= 2 * b
        nc = len(comps)
        for x in range(nc):
            for y in range(x + 1, nc):
                Dx, bx = comps[x]
                Dy, by = comps[y]
                if Dx is Fv and Dy is Fv:
                    continue  # distinct generic fibres are disjoint
                e += bx * by * dot(Dx, Dy)
        # snc / sanity checks on the final configuration
        for i in range(m):
            for j in range(i + 1, m):
                x = dot(Cs[i], Cs[j])
                assert x == self.unblown_pairs.get((i, j), 0), ("strict C.C", i, j, x)
            for k, bl in enumerate(self.blow):
                x = dot(Cs[i], Es[k])
                assert x in (0, 1)
        pos = [k for k, bl in enumerate(self.blow) if bl['b'] > 0]
        for k in pos:
            for l in pos:
                if k < l:
                    assert dot(Es[k], Es[l]) in (0, 1)
        # nefness information
        nefmin = min([dot(KB, D) for D in Cs] + [dot(KB, D) for D in Es] + [Fr(10**9)])
        return dict(KB2=KB2, eorb=e, Q=3 * e - KB2, nefmin=nefmin)

    def ledger(self, S):
        Q = Fr(8, 9) * S + Fr(4, 9) * self.unblown
        Q += sum(4 - Fr(8, 9) * k for k in self.big)
        Q += sum(Fr(12 - 4 * e, 9) for e in self.chains)
        return Q

    def bound(self, S):
        return Fr(8, 9) * S + Fr(4, 9) * self.I - Fr(4, 9) * sum(phi(k) for k in self.big)


def chern_side(m, finite, vinf, degs):
    """finite: list of (weight, v) for finite points; vinf: contact matrix at
    infinity on P^1 x P^1; degs: d_i (C_i^2 = 2 d_i on P^1 x P^1)."""
    for w, v in finite:
        check_ultrametric(v, m)
    check_ultrametric(vinf, m)
    vn, c2, steps = normalise_infinity(m, vinf, [2 * d for d in degs])
    # total intersection numbers on the normalised surface
    Cdot = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i != j:
                Cdot[i][j] = sum(w * v[i][j] for w, v in finite) + vn[i][j]
    for i in range(m):
        for j in range(m):
            if i != j:
                assert 2 * Cdot[i][j] == c2[i] + c2[j], ("2CiCj != Ci^2+Cj^2", i, j)
    S = sum(c2)
    X = Surface(m)
    for w, v in finite:
        X.process_point(v, w)
    X.process_point(vn, 1)
    assert 2 * X.I == (m - 1) * S
    res = X.direct(c2, Cdot)
    res.update(S=S, I=X.I, ledger=X.ledger(S), bound=X.bound(S), steps=steps,
               nblow=len(X.blow), big=sorted(X.big), chains=sorted(X.chains),
               unblown=X.unblown)
    assert res['Q'] == res['ledger'], (res['Q'], res['ledger'])
    assert res['ledger'] <= res['bound']
    return res


def iv2(m, L, layers, alpha, degD, E, delta=None):
    """layers: list of (weight, s_h); alpha: list; degD[(i,j)]; E[(i,j)]."""
    need = max(Fr(degD[(i, j)]) - Fr(alpha[i] + alpha[j], 2) + Fr(L, 4)
               for i in range(m) for j in range(i + 1, m))
    dl = max(Fr(0), need) if delta is None else Fr(delta)
    lhs = Fr(3, 4) * (m - 24) * L + sum(w * (Fr(s) - Fr(m, 2)) ** 2 for w, s in layers)
    rhs = m * (m + 3) * dl + 4 * sum(E.values())
    return dict(L=L, delta=dl, sumE=sum(E.values()), lhs=lhs, rhs=rhs, slack=rhs - lhs)
