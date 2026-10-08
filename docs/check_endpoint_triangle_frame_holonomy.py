"""Exact independent checks for the four-row triangle-frame identities.

The Gaussian fixtures use all 15 nonempty Boolean blocks, with distinct
split-prime norms and corrections that overlap core factors. No thinness
assumption is used. The final fixtures check the independent Pell cycle.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, prod
from typing import Union
import random


Gaussian = tuple[int, int]
Vector = tuple[Union[int, Fraction], Union[int, Fraction]]
Matrix = tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]
LABELS = (1, 2, 3, 4)
FULL_MASK = 15


def mul(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def conj(z: Gaussian) -> Gaussian:
    return z[0], -z[1]


def norm(z: Gaussian) -> int:
    return z[0] * z[0] + z[1] * z[1]


def det(z: Vector, w: Vector) -> int:
    return z[0] * w[1] - z[1] * w[0]


def gprod(values) -> Gaussian:
    out = (1, 0)
    for value in values:
        out = mul(out, value)
    return out


def bit(label: int) -> int:
    return 1 << (label - 1)


def pairmask(i: int, j: int) -> int:
    return bit(i) | bit(j)


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    for d in range(2, isqrt(p) + 1):
        if p % d == 0:
            return False
    return True


def gaussian_split_prime_generators(count: int) -> tuple[Gaussian, ...]:
    out = []
    p = 5
    while len(out) < count:
        if p % 4 == 1 and is_prime(p):
            for a in range(1, isqrt(p) + 1):
                b2 = p - a * a
                b = isqrt(b2)
                if b > 0 and b * b == b2:
                    out.append((a, b))
                    break
            else:
                raise AssertionError(f"split prime {p} has no generator")
        p += 1
    assert len({norm(z) for z in out}) == count
    return tuple(out)


def rows_and_frames(
    blocks: dict[int, Gaussian], corrections: dict[int, Gaussian]
) -> tuple[dict[int, Gaussian], dict[int, Gaussian], dict[int, dict[int, Gaussian]]]:
    rows = {
        i: mul(corrections[i], gprod(blocks[mask] for mask in range(1, 16)
                                     if mask & bit(i)))
        for i in LABELS
    }
    vertices = {
        i: gprod((corrections[i], blocks[bit(i)],
                  conj(blocks[FULL_MASK ^ bit(i)])))
        for i in LABELS
    }
    frames = {}
    for omitted in LABELS:
        triangle = tuple(i for i in LABELS if i != omitted)
        frames[omitted] = {}
        for i in triangle:
            j, k = (v for v in triangle if v != i)
            frames[omitted][i] = gprod((
                vertices[i], blocks[pairmask(i, omitted)],
                conj(blocks[pairmask(j, k)]),
            ))
    return rows, vertices, frames


def common_norm(blocks: dict[int, Gaussian], i: int, j: int) -> int:
    return prod(norm(blocks[mask]) for mask in range(1, 16)
                if mask & bit(i) and mask & bit(j))


def exact_divide(z: Gaussian, n: int) -> Gaussian:
    assert n > 0 and z[0] % n == 0 and z[1] % n == 0
    return z[0] // n, z[1] // n


def matmul(a: Matrix, b: Matrix) -> Matrix:
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))  # type: ignore[return-value]


def matvec(a: Matrix, v: Vector) -> tuple[Fraction, Fraction]:
    return (a[0][0] * v[0] + a[0][1] * v[1],
            a[1][0] * v[0] + a[1][1] * v[1])


def matdet(a: Matrix) -> Fraction:
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def trace(a: Matrix) -> Fraction:
    return a[0][0] + a[1][1]


def columns(v: Vector, w: Vector) -> Matrix:
    return ((Fraction(v[0]), Fraction(w[0])),
            (Fraction(v[1]), Fraction(w[1])))


def inverse(a: Matrix) -> Matrix:
    d = matdet(a)
    assert d
    return ((a[1][1] / d, -a[0][1] / d),
            (-a[1][0] / d, a[0][0] / d))


def transition(frames: dict[int, dict[int, Gaussian]],
               dest: int, source: int) -> Matrix:
    """T_dest,source maps the omitted-source frame to omitted-dest."""
    assert dest != source
    i, j = sorted(set(frames[dest]) & set(frames[source]))
    source_basis = columns(frames[source][i], frames[source][j])
    dest_basis = columns(frames[dest][i], frames[dest][j])
    assert matdet(source_basis) and matdet(dest_basis)
    answer = matmul(dest_basis, inverse(source_basis))
    assert matdet(answer) == 1
    assert matvec(answer, frames[source][i]) == frames[dest][i]
    assert matvec(answer, frames[source][j]) == frames[dest][j]
    return answer


def check_fixture(blocks: dict[int, Gaussian],
                  corrections: dict[int, Gaussian]) -> tuple[bool, int]:
    assert set(blocks) == set(range(1, 16))
    assert set(corrections) == set(LABELS)
    block_norms = tuple(norm(blocks[mask]) for mask in range(1, 16))
    assert all(n > 1 for n in block_norms)
    assert all(gcd(a, b) == 1 for a, b in combinations(block_norms, 2))
    assert all(gcd(abs(z[0]), abs(z[1])) == 1 and (z[0] - z[1]) % 2
               for z in blocks.values())

    rows, vertices, frames = rows_and_frames(blocks, corrections)
    pair_t = {}
    pair_z = {}
    for i, j in combinations(LABELS, 2):
        g = common_norm(blocks, i, j)
        d = det(rows[i], rows[j])
        assert d % g == 0
        pair_t[i, j] = d // g
        pair_z[i, j] = exact_divide(mul(rows[i], conj(rows[j])), g)
        assert pair_z[i, j][1] == -pair_t[i, j]

    if any(t == 0 for t in pair_t.values()):
        return False, 0

    for omitted in LABELS:
        triangle = tuple(i for i in LABELS if i != omitted)
        for i, j in combinations(triangle, 2):
            assert mul(frames[omitted][i], conj(frames[omitted][j])) == pair_z[i, j]
            assert mul(conj(frames[omitted][i]), frames[omitted][j]) == conj(pair_z[i, j])
            assert det(frames[omitted][i], frames[omitted][j]) == pair_t[i, j]
        # The vertex factor is exactly the singleton/triple part of L_i.
        for i in triangle:
            j, k = (v for v in triangle if v != i)
            assert frames[omitted][i] == gprod((
                vertices[i], blocks[pairmask(i, omitted)],
                conj(blocks[pairmask(j, k)]),
            ))

    # The norm equations and additive gcd retain every correction factor.
    a_norms = {i: norm(vertices[i]) for i in LABELS}
    b_norms = (
        norm(blocks[pairmask(1, 2)]) * norm(blocks[pairmask(3, 4)]),
        norm(blocks[pairmask(1, 3)]) * norm(blocks[pairmask(2, 4)]),
        norm(blocks[pairmask(1, 4)]) * norm(blocks[pairmask(2, 3)]),
    )
    b_product = prod(b_norms)
    edge_matching = {(1, 2): 0, (3, 4): 0, (1, 3): 1,
                     (2, 4): 1, (1, 4): 2, (2, 3): 2}

    def oriented_z(i, j):
        return conj(pair_z[i, j]) if i < j else pair_z[j, i]

    def oriented_t(i, j):
        return pair_t[i, j] if i < j else -pair_t[j, i]

    for i, j in combinations(LABELS, 2):
        assert norm(oriented_z(i, j)) == a_norms[i] * a_norms[j] * b_product // b_norms[edge_matching[i, j]]
    for i, j, k in combinations(LABELS, 3):
        assert gprod((oriented_z(i, j), oriented_z(j, k), oriented_z(k, i))) == (
            a_norms[i] * a_norms[j] * a_norms[k] * b_product, 0)
    for i in LABELS:
        eliminants = []
        for j, k in combinations((v for v in LABELS if v != i), 2):
            omitted, = set(LABELS) - {i, j, k}
            zij, zik = oriented_z(i, j), oriented_z(i, k)
            tij, tik, tjk = oriented_t(i, j), oriented_t(i, k), oriented_t(j, k)
            numerator = tuple(tik * zij[d] - tij * zik[d] for d in range(2))
            assert numerator[1] == 0 and numerator[0] % tjk == 0
            quotient = numerator[0] // tjk
            expected = a_norms[i] * norm(blocks[pairmask(i, omitted)]) * norm(blocks[pairmask(j, k)])
            assert quotient == expected == norm(frames[omitted][i])
            eliminants.append(quotient)
        assert sorted(eliminants) == sorted(a_norms[i] * b for b in b_norms)
        assert gcd(gcd(eliminants[0], eliminants[1]), eliminants[2]) == a_norms[i]

    assert (b_norms[0] * pair_t[1, 2] * pair_t[3, 4]
            - b_norms[1] * pair_t[1, 3] * pair_t[2, 4]
            + b_norms[2] * pair_t[1, 4] * pair_t[2, 3]) == 0

    transitions = {(dest, source): transition(frames, dest, source)
                   for dest in LABELS for source in LABELS if dest != source}
    identity = ((Fraction(1), Fraction(0)),
                (Fraction(0), Fraction(1)))
    for dest, source in combinations(LABELS, 2):
        assert matmul(transitions[dest, source], transitions[source, dest]) == identity

    # The order is 4 -> 3 -> 2 -> 4, with T_ab mapping b to a.
    holonomy = matmul(transitions[4, 2],
                      matmul(transitions[2, 3], transitions[3, 4]))
    t = pair_t
    pf = t[1, 2] * t[3, 4] - t[1, 3] * t[2, 4] + t[1, 4] * t[2, 3]
    v1 = frames[4][1]
    v2 = frames[4][2]
    shear = Fraction(-pf, t[1, 3] * t[1, 4])
    expected = matmul(columns(v1, v2), matmul(
        ((Fraction(1), shear), (Fraction(0), Fraction(1))),
        inverse(columns(v1, v2))))
    assert holonomy == expected
    assert matvec(holonomy, v1) == v1
    assert matvec(holonomy, v2) == (
        Fraction(v2[0]) + shear * v1[0],
        Fraction(v2[1]) + shear * v1[1],
    )
    assert matdet(holonomy) == 1
    assert trace(holonomy) == 2

    # A second loop based at frame 4 fixes v2. Check its sign and the
    # product/commutator traces in the same, highly oblique actual basis.
    second = matmul(transitions[4, 1],
                    matmul(transitions[1, 3], transitions[3, 4]))
    basis = columns(v1, v2)
    b_shear = Fraction(pf, t[2, 3] * t[2, 4])
    assert matmul(inverse(basis), matmul(second, basis)) == (
        (Fraction(1), Fraction(0)), (b_shear, Fraction(1)))
    assert trace(matmul(holonomy, second)) == 2 + shear * b_shear
    commutator = matmul(matmul(holonomy, second),
                        matmul(inverse(holonomy), inverse(second)))
    assert trace(commutator) == 2 + (shear * b_shear) ** 2

    # Canonical frames depend on residuals alone. Test the explicit gauge
    # change, denominator budget, and closed-word trace height bound.
    canonical = {}
    gauges = {}
    for omitted in LABELS:
        a, b, c = sorted(frames[omitted])
        canonical[omitted] = {
            a: (Fraction(1), Fraction(0)),
            b: (Fraction(0), Fraction(t[a, b])),
            c: (Fraction(-t[b, c], t[a, b]), Fraction(t[a, c])),
        }
        gauges[omitted] = matmul(
            columns(canonical[omitted][a], canonical[omitted][b]),
            inverse(columns(frames[omitted][a], frames[omitted][b])))
        assert matdet(gauges[omitted]) == 1
        for label in (a, b, c):
            assert matvec(gauges[omitted], frames[omitted][label]) == canonical[omitted][label]
    T = max(map(abs, t.values()))
    L = prod(map(abs, t.values()))
    D = L ** 3
    canonical_transitions = {}
    for dest, source in transitions:
        canonical_transitions[dest, source] = transition(canonical, dest, source)
        direct = canonical_transitions[dest, source]
        assert direct == matmul(gauges[dest], matmul(
            transitions[dest, source], inverse(gauges[source])))
        for row in direct:
            for value in row:
                assert abs(value) <= 2 * T ** 2
                assert D % value.denominator == 0
    for word in ((4, 3, 2, 4), (4, 3, 1, 4), (4, 1, 2, 3, 4)):
        actual_product = canonical_product = identity
        for source, dest in zip(word, word[1:]):
            actual_product = matmul(transitions[dest, source], actual_product)
            canonical_product = matmul(canonical_transitions[dest, source], canonical_product)
        assert trace(actual_product) == trace(canonical_product)
        length = len(word) - 1
        numerator = trace(actual_product) * D ** length
        assert numerator.denominator == 1
        assert abs(numerator) <= 4 ** length * T ** (20 * length)
    return True, pf


def check_boolean_frames() -> tuple[int, int, int]:
    generators = gaussian_split_prime_generators(15)
    rng = random.Random(20261006)
    passed = nonzero_pf = zero_pf = 0
    for trial in range(300):
        shuffled = list(generators)
        rng.shuffle(shuffled)
        blocks = {mask: conj(z) if rng.randrange(2) else z
                  for mask, z in enumerate(shuffled, 1)}
        corrections = {}
        for i in LABELS:
            # The first factor also occurs in P_i; the other chosen factors
            # may overlap either orientation or another row's core.
            values = [blocks[bit(i)]]
            for _ in range(trial % 3 + 1):
                z = blocks[rng.randint(1, 15)]
                values.append(conj(z) if rng.randrange(2) else z)
            corrections[i] = gprod(values)
        ok, pf = check_fixture(blocks, corrections)
        if not ok:
            continue
        passed += 1
        nonzero_pf += pf != 0
        zero_pf += pf == 0
    assert passed >= 250 and nonzero_pf > 0
    return passed, nonzero_pf, zero_pf


def pell_step(u: int, b: int) -> tuple[int, int]:
    return 649 * u + 2340 * b, 180 * u + 649 * b


def check_pell_cycles() -> int:
    u, b = 1, 0
    for n in range(1, 33):
        u, b = pell_step(u, b)
        assert u * u - 13 * b * b == 1
        assert gcd(u, b) == 1 and b % 3 == 0 and u % 3 != 0
        a = (u + 3 * b, 2 * b)
        z = (2 * b, u - 3 * b)
        c = (2 * u + 8 * b, u + b)
        assert c == (2 * a[0] + z[0], 2 * a[1] + z[1])
        assert det(a, z) == 1
        for point in (a, z, c):
            assert gcd(abs(point[0]), abs(point[1])) == 1
            assert (point[0] - point[1]) % 2
        na, nb, nc = map(norm, (a, z, c))
        assert (na, nb, nc) == (
            1 + 26 * b * b + 6 * u * b,
            1 + 26 * b * b - 6 * u * b,
            5 + 130 * b * b + 34 * u * b,
        )
        assert min(na, nb, nc) > 1
        assert all(gcd(x, y) == 1 for x, y in combinations((na, nb, nc), 2))
        zab, zbc, zca = mul(a, conj(z)), mul(z, conj(c)), mul(c, conj(a))
        assert zab == (4 * u * b, -1)
        assert zbc == (1 + 26 * b * b + 2 * u * b, 2)
        assert zca == (2 + 52 * b * b + 16 * u * b, 1)
        assert (3 * zab[0] + 2 * zbc[0] - zca[0],
                3 * zab[1] + 2 * zbc[1] - zca[1]) == (0, 0)
    return 32


def check_isotropic_counterfixture() -> None:
    blocks = dict(enumerate(gaussian_split_prime_generators(15), 1))
    corrections = {i: (1, 0) for i in LABELS}
    ok, pf = check_fixture(blocks, corrections)
    assert ok and pf == 5298313940220
    rows, _, frames = rows_and_frames(blocks, corrections)
    residuals = tuple(det(rows[i], rows[j]) // common_norm(blocks, i, j)
                      for i, j in combinations(LABELS, 2))
    assert residuals == (-2675448, -1081583, -2556377, 3637775, 1624445, -4799520)
    holonomy = matmul(transition(frames, 4, 2),
                      matmul(transition(frames, 2, 3), transition(frames, 3, 4)))
    p = norm(blocks[pairmask(1, 2)])
    assert p == 17 and pf % p == 5

    def reduce_fraction(value: Fraction) -> int:
        assert value.denominator % p
        return value.numerator * pow(value.denominator, -1, p) % p

    reduced = tuple(tuple(reduce_fraction(value) for value in row) for row in holonomy)
    v1 = tuple(value % p for value in frames[4][1])
    z = tuple(value % p for value in frames[4][3])
    hz = tuple(value % p for value in matvec(reduced, z))
    assert reduced == ((16, 2), (15, 3))
    assert v1 == (6, 6) and z == (2, 9) and hz == (16, 6)
    assert norm(z) % p == 0 and norm(hz) % p == 3 and det(z, hz) % p == 4


def check_directed_cut_classification() -> tuple[int, int]:
    edges = tuple((i, j) for i in LABELS for j in LABELS if i != j)

    def exit_count(selected, mask):
        return sum(bool(mask & bit(i)) and not (mask & bit(j)) for i, j in selected)

    pairs = cycles = 0
    for selected in combinations(edges, 2):
        assert any(exit_count(selected, mask) == 1 for mask in range(16))
        pairs += 1
    for selected in combinations(edges, 3):
        has_two_exit_cut = any(exit_count(selected, mask) == 2 for mask in range(16))
        vertices = {v for edge in selected for v in edge}
        is_cycle = len(vertices) == 3 and all(
            sum(i == v for i, _ in selected) == sum(j == v for _, j in selected) == 1
            for v in vertices)
        assert has_two_exit_cut != is_cycle
        cycles += is_cycle
    assert pairs == 66 and cycles == 8
    return pairs, cycles


if __name__ == "__main__":
    passed, nonzero_pf, zero_pf = check_boolean_frames()
    pell = check_pell_cycles()
    check_isotropic_counterfixture()
    pairs, cycles = check_directed_cut_classification()
    print(f"PASS: {passed} nonsingular actual 15-block Gaussian frame fixtures; "
          f"{nonzero_pf} nonzero and {zero_pf} zero Pfaffian shears. "
          "All six Z_ij factorizations, twelve SL2(Fraction) transitions, "
          "two-loop shears and traces, canonical gauges, closed-word height bounds, "
          "seven-norm equations, and four additive gcd identities agree.")
    print(f"PASS: {pell} Pell-cycle fixtures with primitive, unramified, "
          "pairwise norm-coprime A,B,C and imaginary heights (-1,2,1).")
    print(f"PASS: the p=17 isotropic-line counterfixture; all {pairs} directed edge pairs "
          f"and 220 edge triples, with exactly {cycles} directed-cycle exceptions.")
