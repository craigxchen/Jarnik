"""Bounded exact search for four equal-content vertical cofactors.

For fixed ``d`` we enumerate distinct positive ``t`` and set
``w=(-d,t)``.  The Gaussian lcm ``A`` is computed by Euclid in Z[i].
An admissible *actual* endpoint factorization uses ``B=F*A`` with ``B=2P``
coordinatewise even and every ``B/w`` ordinary-primitive.  The minimal useful
multiplier is ``F=1`` when ``A`` is even, or ``F=1+i`` when ``A`` is odd-odd
and every ``A/w`` has odd norm.  If ``A`` has mixed parity, a multiplier making
it even contains 2 and makes every complement non-primitive.  We report all
other primitive lcm candidates as abstract only.

The search is deliberately finite evidence, not a proof.  It covers
exhaustive boxes and arithmetic progressions, plus the Pell orbit recorded in
``vertical_lcm_pell_counterfamily.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from fractions import Fraction
from math import gcd, sqrt
from random import Random

Gaussian = tuple[int, int]


def add(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] + w[0], z[1] + w[1]


def sub(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] - w[0], z[1] - w[1]


def mul(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def conj(z: Gaussian) -> Gaussian:
    return z[0], -z[1]


def norm(z: Gaussian) -> int:
    return z[0] * z[0] + z[1] * z[1]


def divmod_gaussian(z: Gaussian, w: Gaussian) -> tuple[Gaussian, Gaussian]:
    """Euclidean division, using nearest integer real and imaginary parts."""
    den = norm(w)
    # z/w = z*conj(w)/N(w).  Compute nearest integers without float
    # conversion; Pell examples quickly exceed binary64's exact range.
    def nearest(num: int) -> int:
        q, r = divmod(num, den)  # q is floor even when num is negative
        return q + (2 * r >= den)

    qr = nearest(z[0] * w[0] + z[1] * w[1])
    qi = nearest(z[1] * w[0] - z[0] * w[1])
    q = qr, qi
    return q, sub(z, mul(q, w))


def exact_div(z: Gaussian, w: Gaussian) -> Gaussian:
    q, r = divmod_gaussian(z, w)
    if r != (0, 0):
        raise ArithmeticError(f"nonexact Gaussian division {z}/{w}: {q}, {r}")
    return q


def gcd_gaussian(z: Gaussian, w: Gaussian) -> Gaussian:
    while w != (0, 0):
        _, r = divmod_gaussian(z, w)
        z, w = w, r
    return z


def lcm_gaussian(z: Gaussian, w: Gaussian) -> Gaussian:
    if z == (0, 0) or w == (0, 0):
        return (0, 0)
    return mul(exact_div(z, gcd_gaussian(z, w)), w)


def lcm_many(values: list[Gaussian]) -> Gaussian:
    a = (1, 0)
    for z in values:
        a = lcm_gaussian(a, z)
    return a


def ordinary_content(z: Gaussian) -> int:
    return gcd(abs(z[0]), abs(z[1]))


@dataclass(frozen=True)
class Result:
    d: int
    ts: tuple[int, int, int, int]
    A: Gaussian
    complements: tuple[Gaussian, Gaussian, Gaussian, Gaussian]
    primitive: bool
    multiplier: Gaussian | None
    ratio: float

    @property
    def exact_ratio_squared(self) -> tuple[int, int]:
        """( |A|/(min(t/d)^2) )^2 as a reduced fraction."""
        m = min(self.ts)
        num, den = norm(self.A) * self.d**4, m**4
        q = gcd(num, den)
        return num // q, den // q

    @property
    def actual_ratio_squared(self) -> tuple[int, int]:
        """( |B|/(min(t/d)^2) )^2 for the actual integral center."""
        m = min(self.ts)
        num, den = norm(self.actual_B) * self.d**4, m**4
        q = gcd(num, den)
        return num // q, den // q

    @property
    def actual_B(self) -> Gaussian:
        if self.multiplier is None:
            return self.A
        return mul(self.multiplier, self.A)

    @property
    def actual_complements(self) -> tuple[Gaussian, Gaussian, Gaussian, Gaussian]:
        if self.multiplier is None:
            return self.complements
        return tuple(mul(self.multiplier, v) for v in self.complements)

    @property
    def actual(self) -> bool:
        return self.multiplier is not None


def inspect(d: int, ts: tuple[int, int, int, int]) -> Result:
    ws = [(-d, t) for t in ts]
    A = lcm_many(ws)
    vs = tuple(exact_div(A, w) for w in ws)
    primitive = all(ordinary_content(v) == 1 for v in vs)
    multiplier: Gaussian | None = None
    if primitive and A[0] % 2 == 0 and A[1] % 2 == 0:
        # B=A gives P=B/2 integral and retains primitive complements.
        multiplier = (1, 0)
    elif primitive and A[0] % 2 and A[1] % 2:
        # B=(1+i)A is coordinatewise even.  It retains ordinary
        # primitivity only when every (1+i)(A/w) does.
        pi = (1, 1)
        pi_vs = tuple(mul(pi, v) for v in vs)
        if all(ordinary_content(v) == 1 for v in pi_vs):
            multiplier = pi
    B = A if multiplier is None else mul(multiplier, A)
    ratio = sqrt(norm(A)) * d * d / (min(ts) * min(ts))
    return Result(d, ts, A, vs, primitive, multiplier, ratio)


def better(result: Result, best: list[Result], keep: int = 12) -> None:
    if not result.primitive:
        return
    best.append(result)
    best.sort(key=lambda r: Fraction(*r.exact_ratio_squared))
    del best[keep:]


def exhaustive(d: int, tmax: int, keep: int = 12) -> tuple[list[Result], list[Result], int]:
    actual: list[Result] = []
    abstract: list[Result] = []
    count = 0
    for ts in combinations(range(1, tmax + 1), 4):
        r = inspect(d, ts)
        if not r.primitive:
            continue
        count += 1
        if r.actual:
            better(r, actual, keep)
        else:
            better(r, abstract, keep)
    return actual, abstract, count


def arithmetic_progressions(d: int, tmax: int, keep: int = 12) -> tuple[list[Result], int]:
    best: list[Result] = []
    count = 0
    for a in range(1, tmax + 1):
        for h in range(1, (tmax - a) // 3 + 1):
            ts = (a, a + h, a + 2 * h, a + 3 * h)
            r = inspect(d, ts)
            if r.actual:
                count += 1
                better(r, best, keep)
    return best, count


def random_search(d: int, tmax: int, trials: int, seed: int = 0, keep: int = 12) -> list[Result]:
    rng = Random(seed)
    best: list[Result] = []
    for _ in range(trials):
        ts = tuple(sorted(rng.sample(range(1, tmax + 1), 4)))
        r = inspect(d, ts)
        if r.actual:
            better(r, best, keep)
    return best


def advance(z: Gaussian) -> Gaussian:
    x, y = z
    return 17 * x + 4 * y, 4 * x + y


def pell_orbit(steps: int = 20) -> list[tuple[int, tuple[int, int, int, int], Result]]:
    """Four-factor orbit from the companion note, retaining positive t only."""
    a, b, c, e = (0, 1), (1, 3), (4, -3), (2, 1)
    g = (1, 2)
    out = []
    for n in range(steps):
        products = [mul(a, b), mul(a, c), mul(b, c), mul(a, e)]
        ws = [mul((2, 0), mul(g, p)) for p in products]
        ts = tuple(w[1] for w in ws)
        if len(set(ts)) == 4 and min(ts) > 0:
            r = inspect(10, tuple(sorted(ts)))
            # Sorting only changes labels, and lcm/complements are symmetric.
            out.append((n, tuple(sorted(ts)), r))
        a, b, c, e = map(advance, (a, b, c, e))
    return out


def fmt(r: Result) -> str:
    ratio_label = "actual_ratio" if r.actual else "abstract_ratio"
    return (
        f"d={r.d} t={r.ts} A={r.A} N(A)={norm(r.A)} "
        f"B={r.actual_B} multiplier={r.multiplier} "
        f"lcm_ratio={r.ratio:.12g} lcm_ratio^2={r.exact_ratio_squared} "
        f"{ratio_label}={sqrt(Fraction(*r.actual_ratio_squared)):.12g} "
        f"v=A/w={r.complements} v_actual=B/w={r.actual_complements}"
    )


def verify_actual(r: Result) -> None:
    """Check the endpoint circle identity for an actual result."""
    assert r.actual
    B = r.actual_B
    assert B[0] % 2 == 0 and B[1] % 2 == 0
    P = B[0] // 2, B[1] // 2
    ws = [(-r.d, t) for t in r.ts]
    assert all(mul(v, w) == B for v, w in zip(r.actual_complements, ws))
    assert all(ordinary_content(v) == 1 for v in r.actual_complements)
    points = [add(P, (r.d * v[0], r.d * v[1])) for v in r.actual_complements]
    assert len(set(points)) == 4 and all(norm(q) == norm(P) for q in points)


def run() -> None:
    print("TARGETED EXACT CONFIGURATIONS")
    for d, ts in (
        (2, (6, 8, 16, 42)),
        (5, (15, 25, 65, 155)),
        (10, (30, 40, 80, 210)),
    ):
        r = inspect(d, ts)
        verify_actual(r)
        print("  ACTUAL", fmt(r))

    # This box is exhaustive: C(40,4)=91,390 tuples per d.
    print("EXHAUSTIVE SEARCH")
    for d in (1, 2, 5, 10):
        actual, abstract, primitive_count = exhaustive(d, 40)
        print(f"d={d}: primitive lcm tuples={primitive_count}; retained actual={len(actual)} abstract={len(abstract)}")
        for r in actual[:5]:
            print("  ACTUAL", fmt(r))
        for r in abstract[:2]:
            print("  ODD-A", fmt(r))

    print("ARITHMETIC PROGRESSIONS")
    for d in (1, 2, 5, 10):
        best, count = arithmetic_progressions(d, 100)
        print(f"d={d}: actual AP count={count}")
        for r in best[:5]:
            print("  ACTUAL", fmt(r))

    print("RANDOM SEARCH (t<=1000, 1000 trials per d)")
    for d in (1, 2, 5, 10):
        best = random_search(d, 1000, 1_000, seed=100 + d)
        print(f"d={d}: hits={len(best)}")
        for r in best[:5]:
            print("  ACTUAL", fmt(r))

    print("PELL ORBIT (d=10)")
    for n, ts, r in pell_orbit(20):
        print(f"n={n}: {fmt(r)}")


if __name__ == "__main__":
    run()
