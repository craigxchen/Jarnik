"""Exact local valuation certificate; no external dependencies."""
from itertools import combinations, product
import random


def check_case(a, b, size, e, cancellation=None):
    m = len(a)
    tt = set(range(size))
    kappa = [a[i] + b[i] - (e if i in tt else 0) for i in range(m)]
    if min(kappa) < 0:
        return False
    rho = {}
    for i, j in combinations(range(m), 2):
        alpha, beta = b[i] + a[j], a[i] + b[j]
        delta = min(alpha, beta)
        if alpha == beta and cancellation is not None:
            delta += cancellation[i, j]
        core = e if i in tt and j in tt else 0
        if delta < core:
            return False
        rho[i, j] = delta - core
    at = min(a[:size])
    bt = min(b[:size])
    assert e <= at + bt <= e + min(kappa[:size])
    aa, bb = list(a), list(b)
    if size == m:
        # Multiply by (pi/bar(pi))^bt, preserving every row norm.
        aa = [v + bt for v in aa]
        bb = [v - bt for v in bb]
        assert min(bb) == 0
        cc = [min(x, y) for x, y in zip(aa, bb)]
        assert all(cc[i] <= kappa[i] for i in range(m))
        retained = max(0, e - sum(kappa))
        assert min(aa[i] - cc[i] for i in range(m)) >= retained
        assert e - retained <= sum(kappa)
    else:
        cc = [min(x, y) for x, y in zip(aa, bb)]
        cross = sum(v for (i, j), v in rho.items()
                    if (i in tt) != (j in tt))
        assert sum(cc[:size]) <= cross
        assert min(at, bt) <= cross
        assert max(at, bt) >= e - cross
        retained = max(0, e - 2 * cross)
        dominant = aa if at >= bt else bb
        assert min(dominant[i] - cc[i] for i in tt) >= retained
        assert e - retained <= 2 * cross
        assert all(2 * cc[i] <= kappa[i] for i in range(size, m))
    assert sum(cc) <= sum(kappa) + sum(rho.values())
    assert all(min(aa[i] - cc[i], bb[i] - cc[i]) == 0
               for i in range(m))
    return True


def main():
    checked = 0
    for values in product(range(4), repeat=6):
        a, b = values[:3], values[3:]
        for size in range(1, 4):
            for e in range(1, 4):
                checked += check_case(a, b, size, e)
    rng = random.Random(615290)
    for _ in range(30000):
        m = rng.randrange(2, 9)
        a = [rng.randrange(20) for _ in range(m)]
        b = [rng.randrange(20) for _ in range(m)]
        cancellation = {(i, j): rng.randrange(10)
                        for i, j in combinations(range(m), 2)}
        checked += check_case(a, b, rng.randrange(1, m + 1),
                              rng.randrange(1, 15), cancellation)
    # A nonprimitive local pattern repaired by proper-core loss:
    # pi-exponents (4,4,0), conjugate exponents (1,1,0), T={1,2}, e=5.
    assert check_case([4, 4, 0], [1, 1, 0], 2, 5)
    # A full common core split between orientations, removed by one rotation.
    assert check_case([2, 2, 3], [3, 3, 3], 3, 5)
    print(f'Primitive completion verified on {checked} admissible valuation patterns.')


if __name__ == '__main__':
    main()
