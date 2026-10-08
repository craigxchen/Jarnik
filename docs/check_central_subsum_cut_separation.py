"""Exact checks for central subsum separation; no fair-core existence claim."""

from fractions import Fraction
from itertools import combinations
from random import Random

from check_integer_cotangent_affine_content import cleared_clique
from check_integer_cotangent_normalization import primitive_tuple
from check_integer_cotangent_spectral import add, div, mul, neg, power, primitive_vector
from check_six_point_fixed_moment_elliptic import gaussian_valuation, ordinary_valuation


def baseline(es):
    n = len(es)
    h = n // 2
    ordered = sorted(es)
    spread = sum(ordered[h:]) - sum(ordered[:h])
    return [Fraction(sum(abs(e - f) for f in es) - spread, 2) for e in es]


def central(rows):
    nodes = [div(z, rows[0]) for z in rows]
    raw = []
    for i, z in enumerate(nodes):
        denominator = (1, 0)
        for j, w in enumerate(nodes):
            if i != j:
                denominator = mul(denominator, add(z, neg(w)))
        raw.append(div(power(z, len(rows) // 2 - 1), denominator))
    vector, _ = primitive_vector(raw)
    assert all(z[1] == 0 for z in vector) or all(z[0] == 0 for z in vector)
    return [z[0] or z[1] for z in vector]


def small_split_primes():
    answer = []
    for p in range(5, 200, 4):
        if any(p % q == 0 for q in range(2, int(p ** 0.5) + 1)):
            continue
        for a in range(1, int(p ** 0.5) + 1):
            b = int((p - a * a) ** 0.5)
            if b > 0 and a * a + b * b == p:
                answer.append((p, (a, b)))
                break
    return answer


def main():
    rng = Random(20260913)
    prime_data = small_split_primes()
    cliques = valuation_checks = subsum_witnesses = 0
    for n in (4, 6, 8):
        for _ in range(32):
            base = sorted(rng.sample(range(1, 61), n - 1))
            xs, scale = cleared_clique(base)
            rows = primitive_tuple(xs, scale)
            u = central(rows)
            assert all(u) and sum(u) == 0
            for p, pi in prime_data:
                es = [gaussian_valuation(z, pi) for z in rows]
                bs = baseline(es)
                error = (n - 1) * ordinary_valuation(2 * scale, p)
                valuations = [ordinary_valuation(c, p) for c in u]
                assert all(abs(v - b) <= error for v, b in zip(valuations, bs))
                valuation_checks += 1
                levels = sorted(set(es))
                if len(levels) != 2:
                    continue
                for level in levels:
                    minority = [i for i, e in enumerate(es) if e == level]
                    s = len(minority)
                    if s >= n // 2 or (n // 2 - s) * (levels[1] - levels[0]) <= 2 * error:
                        continue
                    for k in range(n):
                        if k in minority:
                            continue
                        assert all(valuations[i] > valuations[k] for i in minority)
                        value = sum(u[i] for i in minority) + u[k]
                        assert value and ordinary_valuation(value, p) == valuations[k]
                        subsum_witnesses += 1
            cliques += 1

    lipschitz_checks = 0
    for n in (4, 6, 8, 10):
        for _ in range(1000):
            es = [rng.randrange(-10, 20) for _ in range(n)]
            fs = [e + rng.randrange(-4, 5) for e in es]
            epsilon = max(abs(e - f) for e, f in zip(es, fs))
            assert all(abs(a - b) <= Fraction(3 * n, 2) * epsilon
                       for a, b in zip(baseline(es), baseline(fs)))
            lipschitz_checks += 1
    assert subsum_witnesses > 0
    print(f'PASS: {cliques} exact circle cliques; {valuation_checks} all-prime central valuation checks.')
    print(f'PASS: {subsum_witnesses} unique-minimum subsum witnesses; {lipschitz_checks} Lipschitz checks.')


if __name__ == '__main__':
    main()
