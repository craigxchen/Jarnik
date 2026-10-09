"""Genuine factor of the C5 determinant on a coordinate 2-plane (x = var_a, y = var_b), others fixed.
Bivariate polys as dict {(i,j): Fraction}. Exact."""
import sys, random, itertools
sys.path.insert(0, '..')
from fractions import Fraction as Fr
from polylib import *
from c5ansatz import LABELS, EDGES, TRIPLES, nbr_edges, triples_of
from c5deg import cleared_matrix

NAMES = ['r%d%d' % e for e in EDGES] + ['s' + ''.join(map(str, T)) for T in TRIPLES]

def params(base, a, b, x, y):
    v = list(base); v[a] = x; v[b] = y
    rho = {e: v[k] for k, e in enumerate(EDGES)}
    sigma = {T: v[5 + k] for k, T in enumerate(TRIPLES)}
    return rho, sigma

def interp1(xs, ys):
    n = len(xs); res = []
    for i in range(n):
        num = [Fr(1)]; den = Fr(1)
        for j in range(n):
            if j != i:
                num = mul(num, linpoly(1, -xs[j])); den *= xs[i] - xs[j]
        res = add(res, scal(ys[i] / den, num))
    return res

def bivariate(base, a, b, N=26):
    # D(x,y) = sum_j c_j(x) y^j ; interpolate in y for each x node, then in x.
    xs = [Fr(3 * k + 1, 2) for k in range(N)]
    ys = [Fr(5 * k - 7, 3) for k in range(N)]
    table = [[det_frac(cleared_matrix(*params(base, a, b, x, y))) for y in ys] for x in xs]
    coeff_in_y = [interp1(ys, row) for row in table]  # for each x: poly in y
    maxdy = max(len(c) for c in coeff_in_y)
    P = {}
    for j in range(maxdy):
        col = [c[j] if j < len(c) else Fr(0) for c in coeff_in_y]
        cx = interp1(xs, col)
        for i, c in enumerate(cx):
            if c != 0:
                P[(i, j)] = c
    return P

def degs(P):
    return max(i for i, j in P), max(j for i, j in P)

def eval2(P, x, y):
    return sum(c * x**i * y**j for (i, j), c in P.items())

def as_poly_in_x(P, y):
    d = {}
    for (i, j), c in P.items():
        d[i] = d.get(i, Fr(0)) + c * y**j
    m = max(d) if d else -1
    return norm([d.get(i, Fr(0)) for i in range(m + 1)])

def divide_linear(P, lin):
    """lin = (alpha, beta, gamma) meaning alpha*x + beta*y + gamma; exact division if possible, else None."""
    al, be, ga = lin
    # divide by treating as poly in x (if al != 0) with coefficients polys in y
    dx, dy = degs(P)
    if al != 0:
        # P = sum_i p_i(y) x^i ; quotient Q with Q*(al x + be y + ga) = P
        p = [norm([P.get((i, j), Fr(0)) for j in range(dy + 1)]) for i in range(dx + 1)]
        q = [None] * dx
        rem = [list(pi) for pi in p]
        for i in range(dx, 0, -1):
            qi = scal(1 / al, rem[i])
            q[i - 1] = qi
            # subtract qi*(be y + ga) x^{i-1}
            rem[i - 1] = sub(rem[i - 1], mul(qi, [ga, be]))
            rem[i] = []
        if norm(rem[0]):
            return None
        Q = {}
        for i, qi in enumerate(q):
            for j, c in enumerate(qi):
                if c != 0:
                    Q[(i, j)] = c
        return Q
    else:
        P2 = {(j, i): c for (i, j), c in P.items()}
        Q = divide_linear(P2, (be, 0, ga))
        return None if Q is None else {(j, i): c for (i, j), c in Q.items()}

if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
    a, b = int(sys.argv[2]), int(sys.argv[3])
    base = [Fr(random.randint(-30, 30), random.randint(1, 4)) for _ in range(10)]
    P = bivariate(base, a, b)
    print('vars', NAMES[a], NAMES[b], 'raw bidegree', degs(P))
    # strip factors x - c, y - c with c among base values / rational roots of slices, and x - y - c
    cands = set(base) | {Fr(0)}
    stripped = []
    changed = True
    while changed:
        changed = False
        # rational roots in x of P(x, y0) common for two y0
        r1 = set(rational_roots(as_poly_in_x(P, Fr(101, 7))))
        r2 = set(rational_roots(as_poly_in_x(P, Fr(-53, 11))))
        for c in r1 & r2:
            Q = divide_linear(P, (Fr(1), Fr(0), -c))
            if Q is not None:
                P = Q; stripped.append(('x', c)); changed = True; break
        if changed: continue
        Py = {(j, i): c for (i, j), c in P.items()}
        r1 = set(rational_roots(as_poly_in_x(Py, Fr(101, 7))))
        r2 = set(rational_roots(as_poly_in_x(Py, Fr(-53, 11))))
        for c in r1 & r2:
            Q = divide_linear(P, (Fr(0), Fr(1), -c))
            if Q is not None:
                P = Q; stripped.append(('y', c)); changed = True; break
        if changed: continue
        for c in [Fr(0)] + [u - v for u in cands for v in cands]:
            for sgn in (1, -1):
                Q = divide_linear(P, (Fr(1), Fr(-sgn), -c))
                if Q is not None:
                    P = Q; stripped.append(('x-%dy' % sgn, c)); changed = True; break
            if changed: break
    print('stripped', len(stripped), [(s, str(c)) for s, c in stripped])
    print('genuine bidegree', degs(P), 'nterms', len(P))
    import pickle
    pickle.dump((base, a, b, P), open('slice_%s_%s.pkl' % (NAMES[a], NAMES[b]), 'wb'))
