"""Exact canonical Gram quotient and source-product checks (stdlib)."""

from itertools import combinations
from random import Random

from check_coupled_orientation_frame_packing import frame, gram
from check_oriented_contact_lattice_rigidity import (
    add, sub, conj, exact, lin, matmul, mul, norm, power, upper,
)


def transpose(M):
    return tuple(zip(*M))


def psi(x, G):
    a, b, c = G[0][0], G[0][1], G[1][1]
    z, w = x
    return add(sub(lin(a, mul(w, w), c, mul(z, z)),
                   lin(2*b, mul(z, w), 0, (0, 0))), (0, 0))


def divisible(z, g):
    z1, n = mul(z, conj(g)), norm(g)
    return z1[0] % n == z1[1] % n == 0


def evaluate(G, lam):
    r, s = lam
    return G[0][0]*r*r+2*G[0][1]*r*s+G[1][1]*s*s


def check_form(U, G, contacts):
    x = frame(U)
    A, B, C = gram(x)
    a, b, c = G[0][0], G[0][1], G[1][1]
    T, xi = C*a-2*B*b+A*c, psi(x, G)
    assert T*T-norm(xi) == 4*(a*c-b*b)
    assert xi[1] % 2 == 0 and (T-xi[0]) % 2 == 0
    H = (((T-xi[0])//2, -xi[1]//2),
         (-xi[1]//2, (T+xi[0])//2))
    assert matmul(matmul(transpose(U), H), U) == G
    delta = (1, 0)
    for generator, lam in contacts:
        exact(lin(lam[0], x[0], lam[1], x[1]), generator)
        assert (evaluate(G, lam) % norm(generator) == 0) == divisible(xi, generator)
        delta = mul(delta, generator)
    in_gamma = all(evaluate(G, lam) % norm(generator) == 0
                   for generator, lam in contacts)
    assert in_gamma == divisible(xi, delta)
    assert norm(xi) <= (a*a+2*b*b+c*c)*(A+C)**2
    if in_gamma and xi != (0, 0):
        assert norm(delta) <= norm(xi)
    return delta


def exact_quotient_and_products():
    U, rng = upper(19), Random(598211)
    x = frame(U)
    specs = [((2, 1), 2), ((3, 2), 2), ((4, 1), 1)]
    contacts = []
    for pi, e in specs:
        value = power(pi, e)
        lam = (value[0]-19*value[1], value[1])
        contacts.append((value, lam))
    delta = (1, 0)
    for generator, _ in contacts:
        delta = mul(delta, generator)
    for _ in range(200):
        a, b, c = [rng.randrange(-100, 101) for _ in range(3)]
        check_form(U, ((a, b), (b, c)), contacts)
        xi = mul(mul((2, 0), delta), (a, b))
        T = 2*c
        H = (((T-xi[0])//2, -xi[1]//2),
             (-xi[1]//2, (T+xi[0])//2))
        G = matmul(matmul(transpose(U), H), U)
        check_form(U, G, contacts)
        assert psi(x, G) == xi
    F = matmul(transpose(U), U)
    assert psi(x, F) == (0, 0)
    for i, j in combinations(range(3), 2):
        M = 1
        for k, (generator, _) in enumerate(contacts):
            if k not in (i, j):
                M *= norm(generator)
        ri, si = contacts[i][1]
        rj, sj = contacts[j][1]
        a, b, c = 2*M*si*sj, -M*(si*rj+ri*sj), 2*M*ri*rj
        G = ((a, b), (b, c))
        check_form(U, G, contacts)
        Qi, Qj = lin(ri, x[0], si, x[1]), lin(rj, x[0], sj, x[1])
        assert psi(x, G) == lin(2*M, mul(Qi, Qj), 0, (0, 0))
        assert a*c-b*b == -M*M*(ri*sj-si*rj)**2
    print('400 forms: exact quotient, reconstruction/parity, full-depth contact equivalence, and ideal/height inequalities pass.')
    print('All three source-product identities and their indefinite determinants pass.')


def cluster_counts():
    for m in range(3, 11):
        A = 2**(m-3)
        clusters = []
        for mask in range(2, 2**m, 2):
            bits = [(mask >> i) & 1 for i in range(m)]
            median = int(sum(bits[:3]) >= 2)
            clusters.append({i for i, b in enumerate(bits) if b != median})
        assert len(clusters) == 4*A-1
        assert len([S for S in clusters if len(S) == 1]) == m
        shared = [S for S in clusters if len(S) != 1]
        for i in range(m):
            assert sum(i in S for S in clusters) == (A if i < 3 else 2*A)
        for i, j in combinations(range(m), 2):
            outside = int(i >= 3)+int(j >= 3)
            overlap = [0, A//2, A][outside]
            assert sum(i in S and j in S for S in clusters) == overlap
            missing = sum(not S.intersection({i, j}) for S in clusters)
            source_cost = outside*A
            assert 2*missing+source_cost == 4*A-2
            missing_shared = sum(not S.intersection({i, j}) for S in shared)
            assert 2*missing_shared+source_cost == 4*A-2*m+2
    print('All private/shared subset-incidence counts and product-height exponents pass for 3 <= m <= 10.')


if __name__ == '__main__':
    exact_quotient_and_products()
    cluster_counts()
