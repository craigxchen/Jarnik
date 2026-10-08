"""Finite shortest-q=1 search for a clean oriented contact lattice.

This is diagnostic only: it does not prove a shortest-vector theorem for
arbitrary contact data.
"""

from itertools import product

from check_oriented_contact_lattice_rigidity import add, conj, hermitian, mul, norm, q


def vector_add(a, b):
    return tuple(add(x, y) for x, y in zip(a, b))


def vector_scale(a, v):
    return tuple(mul(a, x) for x in v)


def vector_combination(a, b, x, y):
    return vector_add(vector_scale(a, x), vector_scale(b, y))


def basis_from_q_one(x, delta):
    y = tuple(
        ((u[0] + mul(delta, conj(u))[0]) // 2,
         (u[1] + mul(delta, conj(u))[1]) // 2)
        for u in x
    )
    assert all((u[0] + mul(delta, conj(v))[0]) % 2 == 0
               and (u[1] + mul(delta, conj(v))[1]) % 2 == 0
               for u, v in zip(x, x))
    assert hermitian(x, y) == (1, 0)
    return y


def search_q_one(x, y, coefficient_bound):
    hits = []
    values = range(-coefficient_bound, coefficient_bound + 1)
    gaussians = list(product(values, repeat=2))
    for a in gaussians:
        for b in gaussians:
            v = vector_combination(a, b, x, y)
            if q(v) == 1:
                hits.append((norm(v[0]) + norm(v[1]), a, b, v))
    return sorted(hits)


def main():
    # Each lambda=(Re G, Im G) makes x=(1,i) satisfy lambda*x=G.
    generators = [(2, 1), (3, 2), (4, 1), (5, 2)]
    delta = (1, 0)
    for generator in generators:
        delta = mul(delta, generator)
    assert norm(delta) == 32045 and norm(delta) % 4 == 1
    x = ((1, 0), (0, 1))
    y = basis_from_q_one(x, delta)
    assert q(x) == 1
    assert hermitian(y, y) == ((1 - norm(delta)) // 2, 0)

    hits = search_q_one(x, y, coefficient_bound=12)
    assert len(hits) == 4
    assert all(a in ((1, 0), (-1, 0), (0, 1), (0, -1)) and b == (0, 0)
               for _, a, b, _ in hits)
    assert hits[0][0] == 2
    print("Clean four-contact lattice: four q=1 vectors in coefficient box 12, all unit multiples of x.")
    print("Observed minimum squared Euclidean norm = 2; this finite search is diagnostic only.")


if __name__ == "__main__":
    main()
