"""Exact CRT, shared-five-vector, and Hermitian basis checks (stdlib)."""

from itertools import combinations
from math import gcd
from random import Random

from check_five_point_cm_norm_lift import bezout, check_system, conj, det, lin, mul, norm
from check_median_rotation_integer_frame import valuation


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def sub(a, b):
    return a[0] - b[0], a[1] - b[1]


def exact(a, b):
    z, n = mul(a, conj(b)), norm(b)
    assert z[0] % n == z[1] % n == 0
    return z[0] // n, z[1] // n


def power(a, n):
    z = (1, 0)
    for _ in range(n):
        z = mul(z, a)
    return z


def matmul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def upper(t):
    return ((1, t), (0, 1))


def lower(t):
    return ((1, 0), (t, 1))


def crt(values, moduli):
    x, modulus = 0, 1
    for a, n in zip(values, moduli):
        x += modulus * ((a - x) * pow(modulus, -1, n) % n)
        modulus *= n
    return x


def iota(pi, exponent):
    p = norm(pi)
    root = (-pi[0] * pow(pi[1], -1, p)) % p
    modulus = p
    for _ in range(1, exponent):
        c = (-(root * root + 1) // modulus * pow(2 * root, -1, p)) % p
        root += c * modulus
        modulus *= p
    assert (root * root + 1) % modulus == 0
    return root


def lift_sl2(a, n):
    assert (a[0][0] * a[1][1] - a[0][1] * a[1][0]) % n == 1
    t = 0
    while gcd(a[0][0] + t * a[1][0], n) != 1:
        t += 1
    u = (a[0][0] + t * a[1][0]) % n
    ui = pow(u, -1, n)
    v = -a[1][0] * ui % n
    b = (a[0][1] + t * a[1][1]) % n
    k = -b * u % n
    matrices = [upper(-t), lower(-v), upper(-k), upper(u), lower(-ui),
                upper(u), ((0, -1), (1, 0))]
    result = ((1, 0), (0, 1))
    for matrix in matrices:
        result = matmul(result, matrix)
    assert det(result[0], result[1]) == 1
    assert all((result[i][j] - a[i][j]) % n == 0 for i in range(2) for j in range(2))
    return result


def complex_det(x, y):
    return sub(mul(x[0], y[1]), mul(x[1], y[0]))


def hermitian(x, y):
    return mul((0, -1), complex_det(tuple(conj(z) for z in x), y))


def q(x):
    return mul(conj(x[0]), x[1])[1]


def check_basis(U, contacts):
    x = ((U[0][0], U[1][0]), (U[0][1], U[1][1]))
    assert q(x) == 1
    delta = (1, 0)
    for generator, vector in contacts:
        exact(lin(vector[0], x[0], vector[1], x[1]), generator)
        delta = mul(delta, generator)
    if delta[0] % 2 == 0:
        delta = mul((0, 1), delta)
    assert delta[0] % 2 == 1 and delta[1] % 2 == 0
    y = []
    for coordinate in x:
        entry = add(coordinate, mul(delta, conj(coordinate)))
        assert entry[0] % 2 == entry[1] % 2 == 0
        y.append((entry[0] // 2, entry[1] // 2))
    y = tuple(y)
    for generator, vector in contacts:
        exact(lin(vector[0], y[0], vector[1], y[1]), generator)
    assert complex_det(x, y) == mul((0, -1), delta)
    D = norm(delta)
    assert hermitian(x, x) == (2, 0)
    assert hermitian(x, y) == (1, 0)
    assert hermitian(y, y) == ((1 - D) // 2, 0)
    for a, b in [((2, 1), (3, -1)), ((-1, 2), (1, 3)), ((0, 0), (1, 0))]:
        vector = tuple(add(mul(a, x[j]), mul(b, y[j])) for j in range(2))
        expected = norm(a) + mul(conj(a), b)[0] + ((1 - D) // 4) * norm(b)
        assert q(vector) == expected
        exact(complex_det(x, vector), delta)
        assert norm(complex_det(x, vector)) <= sum(map(norm, x)) * sum(map(norm, vector))


def random_crt_cases():
    rng = Random(91713)
    primes = [(2, 1), (3, 2), (4, 1)]
    for trial in range(60):
        contacts, matrices, moduli = [], [], []
        for j, pi in enumerate(primes):
            exponent = 1 + (trial + j) % 2
            modulus = norm(pi) ** exponent
            r, s = 0, 0
            while gcd(r, s) != 1:
                r, s = rng.randrange(-20, 21), rng.randrange(-20, 21)
            u, v = bezout(r, s)
            inverse = ((u, v), (-s, r))
            target = ((-iota(pi, exponent), -1), (1, 0))
            matrices.append(matmul(target, inverse))
            moduli.append(modulus)
            contacts.append((power(pi, exponent), (r, s)))
        matrices.append(((1, 0), (0, 1)))
        moduli.append(2)
        n = 1
        for modulus in moduli:
            n *= modulus
        matrix = tuple(tuple(crt([a[i][j] for a in matrices], moduli)
                             for j in range(2)) for i in range(2))
        U = lift_sl2(matrix, n)
        check_basis(U, contacts)
    print('60 arbitrary-slope CRT lifts, including prime powers and parity, pass the integral Hermitian basis checks.')


def shared_five_vector_case():
    primes = [(2, 1), (3, 2), (4, 1), (5, 2), (6, 1),
              (5, 4), (7, 2), (6, 5), (8, 3), (8, 5)]
    edges = list(combinations(range(5), 2))
    data, row_residues, moduli = [], [[] for _ in range(5)], []
    for j, (edge, pi) in enumerate(zip(edges, primes)):
        exponent = 2 if j == 0 else 1
        p = norm(pi)
        mu = int(sum(a in edge for a in range(3)) >= 2)
        cluster = [a for a in range(5) if int(a in edge) != mu]
        outside = [a for a in range(5) if a not in cluster]
        root = iota(pi, exponent + 1)
        allowed = [a for a in range(1, p) if a != 2 * root % p]
        local = {a: p ** exponent * (a + 1) for a in cluster}
        local.update({a: allowed[t] for t, a in enumerate(outside)})
        for a in range(5):
            row_residues[a].append(local[a])
        moduli.append(p ** (exponent + 1))
        data.append((pi, exponent, root, cluster))
    moduli.append(2)
    vs = [(crt(row + [0], moduli), 1) for row in row_residues]
    assert len(set(vs)) == 5
    shear_values, contacts = [], []
    for pi, exponent, root, cluster in data:
        p, pe = norm(pi), norm(pi) ** exponent
        r = vs[cluster[0]][0]
        c = next(t for t in range(p)
                 if all(((vs[a][0] - r) // pe + t) % p for a in cluster))
        shear_values.append(-root - r + pe * c)
        contacts.append((power(pi, exponent), vs[cluster[0]]))
    shear = crt(shear_values + [0], moduli)
    U = upper(shear)
    qs = [(r + shear, s) for r, s in vs]
    for pi, exponent, _, cluster in data:
        for i, row in enumerate(qs):
            assert valuation(row, pi) == (exponent if i in cluster else 0)
            assert valuation(row, conj(pi)) == 0
    assert all(norm(row) % 2 == 1 for row in qs)
    # Put the actual common five-vector system in the required v_1=(1,0)
    # frame before invoking the existing barycentric checker.
    r1 = vs[0][0]
    basis = ((r1, -1), (1, 0))
    U = matmul(U, basis)
    normalized = [(1, r1 - r) for r, _ in vs]
    z, w = (U[0][0], U[1][0]), (U[0][1], U[1][1])
    check_system(z, w, normalized)
    contacts = [(power(pi, exponent), normalized[cluster[0]])
                for pi, exponent, _, cluster in data]
    check_basis(U, contacts)
    print('One shared five-vector system realizes all ten exact clean median contact assignments, including depth two at 5.')
    print('All cross-ratios, barycentric identities, and the canonical Hermitian basis pass; no short-height claim is tested.')


if __name__ == '__main__':
    random_crt_cases()
    shared_five_vector_case()
