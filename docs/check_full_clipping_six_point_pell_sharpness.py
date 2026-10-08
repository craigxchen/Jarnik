"""Exact finite audit of the six-point full-clipping Pell sharpness family."""

from __future__ import annotations

from collections import Counter
from math import isqrt


def mul(z: tuple[int, int], w: tuple[int, int]) -> tuple[int, int]:
    a, b = z
    c, d = w
    return a * c - b * d, a * d + b * c


def conjugate(z: tuple[int, int]) -> tuple[int, int]:
    return z[0], -z[1]


def norm(z: tuple[int, int]) -> int:
    return z[0] * z[0] + z[1] * z[1]


def phase_point(reference: tuple[int, int], h: tuple[int, int]) -> tuple[int, int]:
    numerator = mul(reference, mul(h, h))
    denominator = norm(h)
    assert numerator[0] % denominator == numerator[1] % denominator == 0
    return numerator[0] // denominator, numerator[1] // denominator


def factors(n: int) -> dict[int, int]:
    result: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            result[d] = result.get(d, 0) + 1
            n //= d
        d = 3 if d == 2 else d + 2
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def gaussian_prime(p: int) -> tuple[int, int]:
    for a in range(1, isqrt(p) + 1):
        b = isqrt(p - a * a)
        if b > 0 and a * a + b * b == p:
            return a, b
    raise AssertionError(f"{p} is not split")


def valuation(z: tuple[int, int], pi: tuple[int, int]) -> int:
    a, b = z
    c, d = pi
    p = norm(pi)
    count = 0
    while a or b:
        re, im = a * c + b * d, b * c - a * d
        if re % p or im % p:
            break
        a, b = re // p, im // p
        count += 1
    return count


def clipped_width(source: tuple[int, int], coefficient: tuple[int, int]) -> int:
    left = max(min(source), min(coefficient))
    right = min(max(source), max(coefficient))
    return max(0, right - left)


def audit(b: int, k: int) -> int:
    assert k * k - 2 * b * b == -1
    assert b % 2 == k % 2 == 1
    n = b * b * (b * b + 4)
    z0, w0 = (-2 * b, b * b), (2 * b, b * b)
    hs = [
        (1, 0),
        (-b, 2),
        (k - 2 * b, 1),
        (-k - 2 * b, 1),
        (b + k, 1),
        (b - k, 1),
    ]
    expected_source = [
        (-2 * b, b * b),
        (2 * b, b * b),
        (k, b * b + 1),
        (-k, b * b + 1),
        (-2 * k, b * b - 2),
        (2 * k, b * b - 2),
    ]
    expected_target = [
        (2 * b, b * b),
        (-2 * b, b * b),
        (-2 * k, b * b - 2),
        (2 * k, b * b - 2),
        (k, b * b + 1),
        (-k, b * b + 1),
    ]
    mapped = [(x + b * y, y) for x, y in hs]
    source = [phase_point(z0, h) for h in hs]
    target = [phase_point(w0, h) for h in mapped]
    assert source == expected_source
    assert target == expected_target
    assert set(source) == set(target) and len(set(source)) == 6
    # Return the output to the original source gauge before calling it
    # a self-map. This matrix has projective order two, unlike the shear.
    same_frame = [(-b*x-(b*b+2)*y, 2*x+b*y) for x, y in hs]
    assert [phase_point(z0, h) for h in same_frame] == target
    twice = [(-b*x-(b*b+2)*y, 2*x+b*y) for x, y in same_frame]
    assert twice == [(-(b*b+4)*x, -(b*b+4)*y) for x, y in hs]
    assert all(norm(z) == n for z in source + target)
    twice_stretches = [2 * norm(mh) // norm(h) for h, mh in zip(hs, mapped)]
    assert all(2 * norm(mh) % norm(h) == 0 for h, mh in zip(hs, mapped))
    assert Counter(twice_stretches) == Counter({1: 2, 2: 2, 4: 2})

    u, v = (2, -b), (0, b)
    qclip = 1
    for p, exponent in factors(n).items():
        assert p % 4 == 1
        pi = gaussian_prime(p)
        barpi = conjugate(pi)
        e0 = valuation(z0, pi)
        source_values = [valuation(z, pi) - e0 for z in source]
        source_interval = min(source_values), max(source_values)
        coefficient_interval = (
            valuation(v, pi) - valuation(u, pi),
            valuation(u, barpi) - valuation(v, barpi),
        )
        width = clipped_width(source_interval, coefficient_interval)
        assert width == exponent, (b, p, exponent, source_interval, coefficient_interval)
        assert min(valuation(z, pi) for z in source) == 0
        assert min(valuation(z, barpi) for z in source) == 0
        qclip *= p**width
    assert qclip == n
    assert 4 * n * n // (qclip * qclip) == 4
    return len(factors(n))


def main() -> None:
    b, k = 1, 1
    profiles = 0
    for _ in range(4):
        profiles += audit(b, k)
        k, b = 3 * k + 4 * b, 2 * k + 3 * b
    print(f"PASS: 4 Pell sharpness cases; {profiles} exact split-prime full-clipping profiles")


if __name__ == "__main__":
    main()
