# Referee's independent library (does NOT import any author code).
# Exact arithmetic: Python Fractions for point data, flint.fmpz_poly for Z[s].
from fractions import Fraction as Fr
from itertools import combinations, permutations
import flint
P = flint.fmpz_poly
ZERO = P([0])
LAB = ['O', 'a1', 'a2', 'b1', 'b2', 'c1', 'c2']
IDX = {n: i for i, n in enumerate(LAB)}

# ---- arrangement of L7.md Sec 3.1 (typed from the text, not from author's code)
m1 = (1, 2, -3); m2 = (2, -1, 5); m3 = (-3, 4, 1)
PL = [(0, 0, 0, 1),
      m1 + (2,), m1 + (-3,),
      m2 + (5,), m2 + (-1,),
      m3 + (7,), m3 + (-4,)]

# ---------------- exact rational linear algebra -----------------
def rref_null(rows, n):
    """nullspace basis of a list of rows (Fractions) with n columns"""
    A = [[Fr(x) for x in r] for r in rows]
    piv = []
    r = 0
    for c in range(n):
        pr = None
        for i in range(r, len(A)):
            if A[i][c] != 0:
                pr = i; break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        iv = 1 / A[r][c]
        A[r] = [x * iv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    free = [c for c in range(n) if c not in piv]
    basis = []
    for fc in free:
        v = [Fr(0)] * n
        v[fc] = Fr(1)
        for i, pc in enumerate(piv):
            v[pc] = -A[i][fc]
        basis.append(v)
    return basis

def intvec(v):
    from math import gcd
    from functools import reduce
    den = 1
    for x in v:
        den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
    w = [int(Fr(x) * den) for x in v]
    g = reduce(gcd, [abs(x) for x in w if x != 0], 0)
    return [x // g for x in w] if g else w

def ev(l, x):
    return sum(Fr(a) * Fr(b) for a, b in zip(l, x))

# ---- event points
def meet3(a, b, c):
    ns = rref_null([PL[a], PL[b], PL[c]], 4)
    assert len(ns) == 1
    return intvec(ns[0])

def quad_monos():
    return [(i, j) for i in range(4) for j in range(i, 4)]

def quad_eval_row(x):
    return [Fr(x[i]) * Fr(x[j]) for (i, j) in quad_monos()]

def web(points):
    """basis of quadrics through points, as symmetric 4x4 integer matrices (2*sym form not used)"""
    ns = rref_null([quad_eval_row(x) for x in points], 10)
    mats = []
    for v in ns:
        v = intvec(v)
        # quadric q(x) = sum_{i<=j} v_ij x_i x_j ; symmetric matrix S with x^T S x = 2 q(x)
        S = [[0] * 4 for _ in range(4)]
        for k, (i, j) in enumerate(quad_monos()):
            if i == j:
                S[i][i] += 2 * v[k]
            else:
                S[i][j] += v[k]; S[j][i] += v[k]
        mats.append(S)
    return mats

# ---------------- polynomial linear algebra over Z[s] -----------------
def pdet(M):
    """Bareiss fraction-free determinant of square matrix of fmpz_poly"""
    n = len(M)
    A = [[M[i][j] for j in range(n)] for i in range(n)]
    sign = 1
    prev = P([1])
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = None
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    sw = i; break
            if sw is None:
                return ZERO
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = A[i][j] * A[k][k] - A[i][k] * A[k][j]
                q = num // prev
                assert q * prev == num
                A[i][j] = q
        prev = A[k][k]
    return A[n - 1][n - 1] * sign

def pkernel(M):
    """kernel vector (signed maximal minors) of an r x (r+1) matrix of fmpz_poly, content removed"""
    r = len(M); n = len(M[0]); assert n == r + 1
    K = []
    for j in range(n):
        sub = [[M[i][c] for c in range(n) if c != j] for i in range(r)]
        d = pdet(sub)
        K.append(d if j % 2 == 0 else -d)
    return pnormalize(K)

def pgcd_list(L):
    g = ZERO
    for x in L:
        if x != 0:
            g = x if g == 0 else g.gcd(x)
    return g

def pnormalize(V):
    g = pgcd_list(V)
    if g == 0:
        return V
    W = [x // g for x in V]
    # fix sign so leading nonzero has positive leading coefficient
    return W

def lin(l, v):
    """l: integer 4-vector, v: list of 4 fmpz_poly"""
    r = ZERO
    for a, b in zip(l, v):
        if a != 0:
            r = r + b * a
    return r

def bil(S, x, y):
    r = ZERO
    for i in range(4):
        for j in range(4):
            if S[i][j] != 0:
                r = r + x[i] * y[j] * S[i][j]
    return r

def matsum(coeffs, mats):
    """sum_i coeffs[i] * mats[i], coeffs fmpz_poly, mats integer 4x4 -> 4x4 of fmpz_poly"""
    S = [[ZERO] * 4 for _ in range(4)]
    for c, Mt in zip(coeffs, mats):
        for i in range(4):
            for j in range(4):
                if Mt[i][j] != 0:
                    S[i][j] = S[i][j] + c * Mt[i][j]
    return S

def bilS(S, x, y):
    """S 4x4 fmpz_poly, x,y lists of fmpz_poly"""
    r = ZERO
    for i in range(4):
        for j in range(4):
            if S[i][j] != 0:
                r = r + S[i][j] * x[i] * y[j]
    return r

# ---------------- line through a configuration (Theorem A inverse) -----------------
def line_of_config(V):
    """V: list of 7 homogeneous vectors (v0,v1) of fmpz_poly. Returns p,q (lists of 4 fmpz_poly)
    with L = [p q] (4x2) the kernel of l_a(L v_a)=0. Unknown order: L[i][j] -> index 2i+j."""
    M = []
    for a in range(7):
        row = []
        for i in range(4):
            for j in range(2):
                row.append(V[a][j] * PL[a][i])
        M.append(row)
    K = pkernel(M)
    p = [K[2 * i] for i in range(4)]
    q = [K[2 * i + 1] for i in range(4)]
    return p, q

def bracket(a, b, p, q):
    return lin(PL[a], p) * lin(PL[b], q) - lin(PL[a], q) * lin(PL[b], p)

def dbracket(a, b, p, q, w1, w2):
    la, lb = PL[a], PL[b]
    return (lin(la, w1) * lin(lb, q) + lin(la, p) * lin(lb, w2)
            - lin(la, w2) * lin(lb, p) - lin(la, q) * lin(lb, w1))

# ---------------- configuration brackets -----------------
def det2(u, v):
    return u[0] * v[1] - u[1] * v[0]

FOURS = list(combinations(range(7), 4))

def plcm(a, b):
    return (a * b) // a.gcd(b)

def comp_count(Br, dBr):
    """Br[(a,b)], dBr[(a,b)] for a<b: brackets and their u-derivatives (fmpz_poly in s).
    Returns comp(f,v) of Lemma 5.3 (sum over t in P^1 of min_k ord_t v_k), and #k with W_k != 0."""
    def br(a, b):
        return Br[(a, b)] if a < b else -Br[(b, a)]
    def dbr(a, b):
        return dBr[(a, b)] if a < b else -dBr[(b, a)]
    items = []
    for (a, b, c, d) in FOURS:
        N = br(a, c) * br(b, d); D = br(a, d) * br(b, c)
        Nd = dbr(a, c) * br(b, d) + br(a, c) * dbr(b, d)
        Dd = dbr(a, d) * br(b, c) + br(a, d) * dbr(b, c)
        W = Nd * D - N * Dd
        if W == 0:
            continue
        assert N != 0 and D != 0
        g = N.gcd(D)
        g2 = g * g
        G2 = W.gcd(g2)
        A = W // G2; Bk = g2 // G2
        inf = 2 * max(N.degree(), D.degree()) - W.degree()
        items.append((A, Bk, inf))
    if not items:
        return None, 0
    Bl = P([1])
    for (A, Bk, inf) in items:
        Bl = plcm(Bl, Bk)
    G = ZERO
    for (A, Bk, inf) in items:
        t = A * (Bl // Bk)
        G = t if G == 0 else G.gcd(t)
    fin = G.degree() - Bl.degree()
    inf = min(it[2] for it in items)
    return fin + inf, len(items)

# ---------------- F_A tangent along a test curve given as configuration V -----------------
E1 = [3, -1, 4, 1]; E2 = [-2, 5, 1, 7]   # fixed generic complement vectors (referee's choice)

def FA_brackets(V, webm):
    p, q = line_of_config(V)
    rows = [[bil(Qm, p, p) for Qm in webm],
            [bil(Qm, q, q) for Qm in webm],
            [bil(Qm, p, q) for Qm in webm]]
    c = pkernel(rows)
    assert any(x != 0 for x in c), "web quadric not unique along curve"
    B = matsum(c, webm)
    e1 = [P([x]) for x in E1]; e2 = [P([x]) for x in E2]
    bp1, bp2 = bilS(B, p, e1), bilS(B, p, e2)
    bq1, bq2 = bilS(B, q, e1), bilS(B, q, e2)
    T = [[bp1, bp2, ZERO, ZERO], [ZERO, ZERO, bq1, bq2], [bq1, bq2, bp1, bp2]]
    k = pkernel(T)
    assert any(x != 0 for x in k)
    w1 = [k[0] * E1[i] + k[1] * E2[i] for i in range(4)]
    w2 = [k[2] * E1[i] + k[3] * E2[i] for i in range(4)]
    Br = {}; dBr = {}
    for a, b in combinations(range(7), 2):
        Br[(a, b)] = bracket(a, b, p, q)
        dBr[(a, b)] = dbracket(a, b, p, q, w1, w2)
    return Br, dBr, (p, q, c, B, k)

def forget_brackets(V, i, w=(1, 0)):
    wv = [P([w[0]]), P([w[1]])]
    Br = {}; dBr = {}
    for a, b in combinations(range(7), 2):
        Br[(a, b)] = det2(V[a], V[b])
        d = ZERO
        if a == i:
            d = d + det2(wv, V[b])
        if b == i:
            d = d + det2(V[a], wv)
        dBr[(a, b)] = d
    return Br, dBr

# ---------------- class of a configuration curve -----------------
def class_of_config(V):
    """Returns dict S(frozenset)->D_S degree, plus flags. Requires 'good'-type genericity checks:
    each point on at most one boundary divisor; returns None if a check fails.
    Works for curves whose collisions are only of pair/triple type at distinct parameters,
    with nothing at s=infinity."""
    br = {}
    for a, b in combinations(range(7), 2):
        br[(a, b)] = det2(V[a], V[b])
    # all brackets degree must be full (no collision at infinity): we check by homogeneous degree
    return br

def class_and_goodness(V):
    """V: 7 homogeneous vectors of fmpz_poly in s. Returns (cls, good, info):
    cls: dict frozenset(S) -> degree (only |S| in {2,3}); good: True iff every root (incl. s=inf)
    has exactly one collision cluster, of size 2 or 3, with simple vanishing of each of its brackets."""
    br = {(a, b): det2(V[a], V[b]) for a, b in combinations(range(7), 2)}
    hd = max(max(V[a][0].degree(), V[a][1].degree()) for a in range(7))  # homogeneous degree of each V_a
    # infinity: bracket has homogeneous degree 2*hd; root at inf of multiplicity 2hd - deg
    for k, v in br.items():
        if v == 0:
            return None, False, 'identically colliding pair %s' % (k,)
    # remove a common factor of all brackets (artefact of a degenerating basis), incl. a common root at infinity
    Gc = pgcd_list(list(br.values()))
    br = {k: v // Gc for k, v in br.items()}
    e = max(v.degree() for v in br.values())
    inf_pairs = [k for k, v in br.items() if v.degree() < e]
    if inf_pairs:
        return None, False, 'collision at infinity %s' % inf_pairs
    prod = P([1])
    for v in br.values():
        prod = plcm(prod, v)
    fac = prod.factor()[1]
    cls = {}
    good = True
    info = []
    for (phi, e) in fac:
        if phi.degree() == 0:
            continue
        pairs = []
        for k, v in br.items():
            m = 0; w = v
            while True:
                qq, rr = divmod(w, phi)
                if rr == 0:
                    m += 1; w = qq
                else:
                    break
            if m:
                pairs.append((k, m))
        labs = set()
        for (k, m) in pairs:
            labs |= set(k)
        # clusters via union find
        parent = {x: x for x in labs}
        def find(x):
            while parent[x] != x:
                x = parent[x]
            return x
        for (k, m) in pairs:
            parent[find(k[0])] = find(k[1])
        clusters = {}
        for x in labs:
            clusters.setdefault(find(x), set()).add(x)
        cl = list(clusters.values())
        mults = [m for (k, m) in pairs]
        if len(cl) != 1 or len(cl[0]) > 3 or max(mults) > 1 or len(pairs) != len(cl[0]) * (len(cl[0]) - 1) // 2:
            good = False
            info.append((str(phi), pairs))
            continue
        S = frozenset(cl[0])
        cls[S] = cls.get(S, 0) + phi.degree()
    return cls, good, info

def class_vector(cls):
    """56-vector over the boundary divisors indexed by subsets S with |S|=2 or 3 (D_S=D_{S^c})"""
    keys = [frozenset(c) for c in combinations(range(7), 2)] + [frozenset(c) for c in combinations(range(7), 3)]
    return tuple(cls.get(k, 0) for k in keys)

# ---------------- member of V_A through a configuration point x -----------------
def member_through(x, webm):
    """x: 7 rationals (affine positions). Returns V(s): 7 homogeneous linear vectors in s (Mobius curve)
    and data. R_x, Q_x, ruling through R_x, parametrised by lines meeting two opposite-ruling lines."""
    V0 = [[P([int(Fr(xa).numerator)]), P([int(Fr(xa).denominator)])] for xa in x]
    p, q = line_of_config(V0)
    p = [int(z.coeffs()[0]) if z != 0 else 0 for z in p]
    q = [int(z.coeffs()[0]) if z != 0 else 0 for z in q]
    # web quadric containing R_x
    def b(S, u, v):
        return sum(Fr(S[i][j]) * u[i] * v[j] for i in range(4) for j in range(4))
    rows = [[b(Qm, p, p) for Qm in webm], [b(Qm, q, q) for Qm in webm], [b(Qm, p, q) for Qm in webm]]
    ns = rref_null(rows, 4)
    assert len(ns) == 1, "quadric through R_x not unique"
    c = intvec(ns[0])
    S = [[sum(c[k] * webm[k][i][j] for k in range(4)) for j in range(4)] for i in range(4)]
    # smoothness
    import sympy
    detS = sympy.Matrix(S).det()
    assert detS != 0, "Q_x singular"
    # opposite-ruling line through y in R: R'_y = span(y, u)
    def opp(y, y2):
        # w0 with B(y,w0)=0, not in span(y,y2)
        cands = rref_null([[sum(S[i][j] * y[j] for j in range(4)) for i in range(4)]], 4)
        for w0 in cands:
            w0 = intvec(w0)
            if sympy.Matrix([y, y2, w0]).rank() == 3:
                break
        Qw0 = b(S, w0, w0) / 2
        By2w0 = b(S, y2, w0) / 2
        assert By2w0 != 0
        u = [Qw0 * y2[i] - 2 * By2w0 * w0[i] for i in range(4)]
        assert b(S, u, u) == 0 and b(S, y, u) == 0
        return intvec(u)
    y0 = p; y1 = q
    u0 = opp(y0, y1); u1 = opp(y1, y0)
    # z(s) = y0 + s u0 ; z'(s) = B(z,u1) y1 - B(z,y1) u1
    zs = [P([y0[i], u0[i]]) for i in range(4)]
    Bzu1 = P([int(b(S, y0, u1)), int(b(S, u0, u1))])
    Bzy1 = P([int(b(S, y0, y1)), int(b(S, u0, y1))])
    zp = [Bzu1 * y1[i] - Bzy1 * u1[i] for i in range(4)]
    V = []
    for a in range(7):
        V.append([lin(PL[a], zp), -lin(PL[a], zs)])
    return V, dict(p=p, q=q, c=c, S=S, u0=u0, u1=u1)
