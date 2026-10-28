"""Referee r3 library (independent of the author's code).
Exact gain of a pencil f = lam * A / C on Mbar_{0,8}, Kapranov chart (anchor 8 at infinity, coords x0..x6).
A, C: arbitrary coprime homogeneous translation-invariant integer polynomials of equal degree (units or not).
All fibre data computed on X from exact factorizations over Z (flint) and exact orders nu_T along L_T
(substitution x_i -> x_{i0} + t x_i, i in T minus i0; lowest power of t).
Output: L.beta computed from F_0, F_1 and F_inf separately (must agree), gain, n(f), N(u) etc.
"""
import itertools
from fractions import Fraction
import flint

NV = 7
CTX = flint.fmpz_mpoly_ctx.get(('x', NV + 1))   # x0..x6, x7 = t
G = CTX.gens()
X = G[:NV]
T_ = G[NV]
PAIRS = list(itertools.combinations(range(NV), 2))
BIG = [T for r in range(3, 7) for T in itertools.combinations(range(NV), r)]   # 98 centres L_T
ALLT = [frozenset(p) for p in PAIRS] + [frozenset(T) for T in BIG]           # 119 divisors
assert len(ALLT) == 119

_subst = {}
def subst_list(T):
    if T not in _subst:
        T = tuple(sorted(T))
        i0 = T[0]
        s = []
        for i in range(NV):
            s.append(X[i0] + T_ * X[i] if (i in T and i != i0) else X[i])
        s.append(T_)
        _subst[T] = s
    return _subst[T]

def tdeg_min(P):
    return min(m[NV] for m in P.monoms())

def nu(P, T):
    """order of homogeneous P along L_T (|T|>=3) or exponent of x_i - x_j (|T|=2)"""
    T = tuple(sorted(T))
    if len(T) == 2:
        i, j = T
        d = X[i] - X[j]
        k = 0
        Q = P
        while True:
            q, r = divmod(Q, d)
            if r != 0:
                return k
            Q = q; k += 1
    return tdeg_min(P.compose(*subst_list(T)))

def initial_form(P, T):
    """lowest-t coefficient of P after substitution, as a polynomial (t set to 1 on that part)"""
    Q = P.compose(*subst_list(T))
    m0 = tdeg_min(Q)
    terms = {}
    for mon, c in zip(Q.monoms(), Q.coeffs()):
        if mon[NV] == m0:
            terms[tuple(mon[:NV]) + (0,)] = c
    return CTX.from_dict(terms)

def diff(i, j):
    return X[i] - X[j]

def is_pair(Y):
    """return frozenset pair if Y = +-(x_i - x_j)"""
    if Y.total_degree() != 1:
        return None
    for (i, j) in PAIRS:
        d = diff(i, j)
        if Y == d or Y == -d:
            return frozenset((i, j))
    return None

_fcache = {}
def factor_data(P):
    """list of (kind, data, mult): kind 'pair' with frozenset, or 'int' with (deg, beta-degree, {T: nu})"""
    c, facs = P.factor()
    out = []
    for Y, m in facs:
        pr = is_pair(Y)
        if pr is not None:
            out.append(('pair', pr, int(m)))
            continue
        key = str(Y)
        if key not in _fcache:
            nus = {frozenset(T): int(nu(Y, T)) for T in BIG}
            deg = int(Y.total_degree())
            _fcache[key] = (deg, 31 * deg - sum(nus.values()), nus)
        out.append(('int', _fcache[key], int(m)))
    return out

def nuvec(fd):
    """nu_T of the product, for all 119 T"""
    v = {T: 0 for T in ALLT}
    for kind, dat, m in fd:
        if kind == 'pair':
            for T in ALLT:
                if dat <= T:
                    v[T] += m
        else:
            for T, n in dat[2].items():
                v[T] += m * n
    return v

def fibre(fd, nuv, mu):
    """beta-degree with multiplicity, reduced beta-degree, interior reduced beta-degree, boundary set"""
    tot = 0; red_int = 0; bd = set()
    for kind, dat, m in fd:
        if kind == 'int':
            tot += m * dat[1]; red_int += dat[1]
    for T in ALLT:
        o = nuv[T] - mu[T]
        assert o >= 0
        if o > 0:
            tot += o; bd.add(T)
    return tot, red_int, bd

def gain_data(A, C, lam):
    """lam: Fraction. f = lam*A/C. B := q*C - p*A  (zero set of f - 1)."""
    lam = Fraction(lam)
    p, q = lam.numerator, lam.denominator
    B = q * C - p * A
    assert B != 0
    fA, fC, fB = factor_data(A), factor_data(C), factor_data(B)
    nA, nC, nB = nuvec(fA), nuvec(fC), nuvec(fB)
    mu = {T: min(nA[T], nC[T]) for T in ALLT}
    for T in ALLT:
        assert nB[T] >= mu[T]
    L0, R0i, bd0 = fibre(fA, nA, mu)
    Li, Rii, bdi = fibre(fC, nC, mu)
    L1, R1i, bd1 = fibre(fB, nB, mu)
    assert L0 == Li == L1, ("fibre classes disagree", L0, Li, L1)
    assert not (bd0 & bd1) and not (bd0 & bdi) and not (bd1 & bdi)
    n = len(bd0) + len(bd1) + len(bdi)
    Fred = R0i + R1i + Rii + n
    gain = L0 - Fred
    # multiplicity form: sum ((m-3)/3) Y.beta must equal gain
    mf = Fraction(0)
    for fd, nv in ((fA, nA), (fC, nC), (fB, nB)):
        for kind, dat, m in fd:
            if kind == 'int':
                mf += Fraction(m - 3, 3) * dat[1]
        for T in ALLT:
            o = nv[T] - mu[T]
            if o > 0:
                mf += Fraction(o - 3, 3)
    assert mf == gain, (mf, gain)
    return dict(L=L0, gain=gain, n=n, n0=len(bd0), n1=len(bd1), ninf=len(bdi),
                R_int=(R0i, R1i, Rii), Bfactors=[(k, (sorted(d) if k == 'pair' else (d[0], d[1])), m) for k, d, m in fB])

def unit_poly(exps):
    """exps: dict pair->exponent >=0 ; returns product"""
    P = CTX.from_dict({(0,) * (NV + 1): 1})
    for (i, j), e in exps.items():
        P = P * diff(i, j) ** e
    return P

def kf_unit(e):
    """k(f) for a unit with exponent vector e (dict pair -> int)"""
    I = set(); deg = [0] * NV
    for (i, j), v in e.items():
        if v:
            I |= {i, j}; deg[i] += v; deg[j] += v
    return len(I) + (1 if any(deg) else 0)

LAMBDA = {4: 16, 5: 24, 6: 28, 7: 30, 8: 31}

def special_lambdas(A, C, rng, nlines=2):
    """candidate lam (Fraction) where B = C - lam A (up to scaling) can be special:
    (i) initial forms along some L_T / pair with equal order proportional;
    (ii) rational roots of disc_s(C(P+sQ) - lam A(P+sQ)) for nlines random lines (union)."""
    cands = set()
    for T in ALLT:
        Tt = tuple(sorted(T))
        if nu(A, Tt) != nu(C, Tt):
            continue
        if len(Tt) == 2:
            i, j = Tt
            sub = [X[k] if k != j else X[i] for k in range(NV)] + [T_]
            a, c = A.compose(*sub), C.compose(*sub)
        else:
            a, c = initial_form(A, Tt), initial_form(C, Tt)
        # proportional?  c = lam a
        ma, mc = dict(zip(map(tuple, a.monoms()), a.coeffs())), dict(zip(map(tuple, c.monoms()), c.coeffs()))
        if set(ma) != set(mc):
            continue
        k0 = next(iter(ma))
        lam = Fraction(int(mc[k0]), int(ma[k0]))
        if all(Fraction(int(mc[k]), int(ma[k])) == lam for k in ma):
            cands.add(lam)
    ctx2 = flint.fmpz_mpoly_ctx.get(('s', 'l'))
    s, l = ctx2.gens()
    good = 0; tries = 0
    while good < nlines:
        tries += 1
        if tries > 50:
            raise RuntimeError('no good line')
        P = [rng.randint(-50, 50) for _ in range(NV)]
        Q = [rng.randint(-50, 50) for _ in range(NV)]
        sub = [P[k] + Q[k] * s for k in range(NV)] + [ctx2.from_dict({(0, 0): 0})]
        # compose into ctx2: build via evaluating monomials
        def restrict(F):
            out = ctx2.from_dict({(0, 0): 0})
            for mon, cf in zip(F.monoms(), F.coeffs()):
                term = ctx2.from_dict({(0, 0): int(cf)})
                for k in range(NV):
                    if mon[k]:
                        term = term * sub[k] ** mon[k]
                out = out + term
            return out
        a, c = restrict(A), restrict(C)
        Gp = c - l * a
        R = Gp.resultant(Gp.derivative(0), 's')   # polynomial in l
        if R == 0 or Gp.degrees()[0] != max(A.total_degree(), C.total_degree()):
            continue   # degenerate line (degree drop in s): draw another
        good += 1
        coeffs = {}
        for mon, cf in zip(R.monoms(), R.coeffs()):
            assert mon[0] == 0
            coeffs[mon[1]] = int(cf)
        up = flint.fmpz_poly([coeffs.get(i, 0) for i in range(max(coeffs) + 1)])
        cst, facs = up.factor()
        for Fct, m in facs:
            if Fct.degree() == 1:
                a0, a1 = int(Fct[0]), int(Fct[1])
                cands.add(Fraction(-a0, a1))
    cands.discard(Fraction(0))
    return cands
