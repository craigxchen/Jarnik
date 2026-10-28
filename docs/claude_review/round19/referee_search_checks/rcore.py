"""Independent referee core (does not import the author's code).
Labels 0..6 = x_1..x_7. Edge (i,j), i<j, polynomial x_i - x_j.
Gain = theory.md Prop 2.6 formula, implemented literally:
  L.beta = 31 D - sum_{|T|>=3} mu_T
  F_red.beta = #{T visible: max(nu_T(A),nu_T(B),nu_T(C)) > mu_T} + sum_{Y interior factor of B} (31 deg Y - sum_{|T|>=3} nu_T(Y))
"""
import itertools, random
from fractions import Fraction
import numpy as np
import flint

N = 7
EDGES = [(i, j) for i in range(N) for j in range(i + 1, N)]
EIDX = {e: k for k, e in enumerate(EDGES)}
VIS = [frozenset(c) for r in range(2, 7) for c in itertools.combinations(range(N), r)]  # 119
BIG = [T for T in VIS if len(T) >= 3]  # 98
assert len(VIS) == 119 and len(BIG) == 98
# incidence: edge inside T
INC = np.array([[1 if (e[0] in T and e[1] in T) else 0 for e in EDGES] for T in VIS], dtype=np.int64)

CTX = flint.fmpz_mpoly_ctx.get(tuple('x%d' % i for i in range(N)) + ('e',), 'lex')
G = CTX.gens()
X = G[:N]
E = G[N]


def vec(edges):
    v = np.zeros(21, dtype=np.int64)
    for e in edges:
        v[EIDX[tuple(e)]] += 1
    return v


def N0(av, cv):
    return int(np.count_nonzero(INC @ (av - cv)))


def poly_of(edges):
    p = 1
    for (i, j) in edges:
        p = p * (X[i] - X[j])
    return p


def nu_T(F, T):
    """exact order of F along L_T (|T|>=3): x_i -> x_i0 + e*x_i, i in T minus i0; lowest e-degree."""
    T = sorted(T)
    i0 = T[0]
    imgs = list(X) + [E]
    for i in T[1:]:
        imgs[i] = X[i0] + E * X[i]
    Fs = F.compose(*imgs)
    return int(min(int(m[N]) for m in Fs.monoms()))


def classify_factor(Y):
    """return ('pair',(i,j),sign) if Y = +-(x_i - x_j), else ('int',)"""
    if Y.total_degree() == 1:
        d = Y.to_dict()
        if len(d) == 2 and sorted(abs(int(c)) for c in d.values()) == [1, 1] and sum(int(c) for c in d.values()) == 0:
            idx = sorted(m.index(1) for m in d.keys())
            return ('pair', tuple(idx))
    return ('int',)


def factorZ(B):
    """flint factorisation; python-flint 0.9 sort-key overflow -> independent sympy fallback."""
    try:
        return B.factor()
    except OverflowError:
        import sympy
        xs = sympy.symbols('x0:7')
        P = sympy.Poly(sympy.sympify(str(B).replace('^', '**')), *xs)
        c, fl = P.factor_list()
        out = []
        for (Fp, m) in fl:
            dd = {tuple(k) + (0,): int(v) for k, v in Fp.as_dict().items()}
            out.append((CTX.from_dict(dd), int(m)))
        # sanity: product reproduces B
        prod = CTX.from_dict({(0,) * 8: int(c)})
        for (Fp, m) in out:
            prod = prod * Fp ** m
        assert prod == B
        return c, out


def full_gain(Aedges, Cedges, p, q, detail=False):
    """lambda = p/q; B = q C - p A.  Literal Prop 2.6 formula, all factors over Z (flint)."""
    D = len(Aedges)
    av, cv = vec(Aedges), vec(Cedges)
    aT = INC @ av
    cT = INC @ cv
    A = poly_of(Aedges)
    C = poly_of(Cedges)
    B = q * C - p * A
    assert not B.is_zero()
    content, facs = factorZ(B)
    nuB = {T: 0 for T in VIS}
    interior = []
    for (Y, m) in facs:
        m = int(m)
        if Y.total_degree() == 0:
            continue
        cl = classify_factor(Y)
        if cl[0] == 'pair':
            nuB[frozenset(cl[1])] += m
            for T in BIG:
                if cl[1][0] in T and cl[1][1] in T:
                    nuB[T] += m
        else:
            nus = {T: nu_T(Y, T) for T in BIG}
            for T in BIG:
                nuB[T] += m * nus[T]
            ybeta = 31 * Y.total_degree() - sum(nus.values())
            interior.append((str(Y), m, Y.total_degree(), ybeta))
    mu = {}
    for k, T in enumerate(VIS):
        mu[T] = min(int(aT[k]), int(cT[k])) if len(T) >= 3 else 0
    Lb = 31 * D - sum(mu[T] for T in BIG)
    cnt = 0
    for k, T in enumerate(VIS):
        if max(int(aT[k]), nuB[T], int(cT[k])) > mu[T]:
            cnt += 1
    Fb = cnt + sum(yb for (_, _, _, yb) in interior)
    g = int(Lb - Fb)
    if detail:
        return g, dict(L=Lb, F=Fb, nbound=cnt, interior=interior, N0=N0(av, cv))
    return g


# ---------- line restriction tools ----------
def lin(p, q, e):
    i, j = e
    return flint.fmpz_poly([p[i] - p[j], q[i] - q[j]])


def restrict(edges, p, q):
    r = flint.fmpz_poly([1])
    for e in edges:
        r = r * lin(p, q, e)
    return r


def sqfree_form(b, D):
    """b(t) = dehomogenised binary form of formal degree D; squarefree on P^1?"""
    if b.is_zero():
        return False
    if b.degree() < D - 1:
        return False
    return b.gcd(b.derivative()).degree() == 0


def disc_poly(a, c, D):
    """Delta(lam) = Disc_D(c - lam a) by exact interpolation at 2D-1 integer lam
    where the leading coefficient does not vanish (then affine disc = formal Disc_D)."""
    pts = []
    lam0 = -D
    aD = a.coeffs()[D] if a.degree() == D else 0
    cD = c.coeffs()[D] if c.degree() == D else 0
    while len(pts) < 2 * D - 1:
        if cD - lam0 * aD != 0:
            f = c - lam0 * a
            pts.append((lam0, int(f.discriminant())))
        lam0 += 1
    # Newton interpolation over Q
    xs = [Fraction(x) for x, _ in pts]
    cf = [Fraction(y) for _, y in pts]
    n = len(xs)
    for k in range(1, n):
        for i in range(n - 1, k - 1, -1):
            cf[i] = (cf[i] - cf[i - 1]) / (xs[i] - xs[i - k])
    poly = [Fraction(0)] * n
    # expand
    res = [Fraction(0)]
    for k in range(n - 1, -1, -1):
        # res = res*(x - xs[k]) + cf[k]
        new = [Fraction(0)] * (len(res) + 1)
        for i, r in enumerate(res):
            new[i + 1] += r
            new[i] -= r * xs[k]
        new[0] += cf[k]
        res = new
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    for r in res:
        assert r.denominator == 1
    return flint.fmpz_poly([int(r) for r in res])


def rational_roots(P):
    out = set()
    if P.is_zero():
        return None
    _, facs = P.factor()
    for (f, m) in facs:
        if f.degree() == 1:
            b, a = int(f.coeffs()[0]), int(f.coeffs()[1])
            fr = Fraction(-b, a)
            out.add(fr)
    return out
