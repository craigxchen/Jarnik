"""Finite LP audit for simultaneous Walsh coset-character heights.

The variables are class weights W[a], a in F_2^t minus {0}.  For every
subspace V and nonzero restriction ell, k is maximized over all row cosets,
and the exact one-flip condition is

    |V| * sum_(a|V=ell) W[a] - sum_a W[a] >= 4*k - 2*|V|.

The scipy solve is only a locator.  The reported optimum is accepted only
after reconstructing a rational primal and dual pair and checking every
inequality exactly.  This file deliberately uses bin(...).count rather
than int.bit_count for compatibility with the older audit environment.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
from math import gcd
from pathlib import Path
import random
import sys


def parity(x: int) -> int:
    return bin(x).count("1") & 1


def dot(a: int, x: int) -> int:
    return parity(a & x)


def subspace_basis(t: int, dimension: int):
    """Enumerate every subspace once as an RREF binary basis tuple."""
    if dimension == 0:
        yield ()
        return
    positions = range(t)
    for pivots in combinations(positions, dimension):
        free = []
        for i, p in enumerate(pivots):
            free.append([q for q in positions if q > p and q not in pivots])
        slots = sum(len(x) for x in free)
        for bits in range(1 << slots):
            basis = []
            offset = 0
            for i, p in enumerate(pivots):
                mask = 1 << p
                for q in free[i]:
                    if (bits >> offset) & 1:
                        mask |= 1 << q
                    offset += 1
                basis.append(mask)
            yield tuple(basis)


def subspace_elements(basis):
    values = [0]
    for v in basis:
        values += [x ^ v for x in values]
    return tuple(sorted(values))


def coset_representatives(elements, m):
    unseen = set(range(m))
    reps = []
    while unseen:
        x = min(unseen)
        coset = {x ^ v for v in elements}
        reps.append(x)
        unseen -= coset
    return reps


def map_random(m, seed):
    labels = list(range(1, m)) + [1]
    random.Random(seed).shuffle(labels)
    return labels


M32 = (
    28, 21, 11, 30, 24, 23, 13, 3,
    12, 7, 25, 4, 31, 1, 22, 3,
    19, 18, 10, 9, 2, 16, 17, 27,
    26, 6, 15, 5, 14, 8, 20, 29,
)


def constraints(labels):
    m = len(labels)
    t = m.bit_length() - 1
    variables = list(range(1, m))
    index = {a: i for i, a in enumerate(variables)}
    rows = []
    metadata = []
    for d in range(1, t + 1):
        for basis in subspace_basis(t, d):
            elements = subspace_elements(basis)
            h = len(elements)
            reps = coset_representatives(elements, m)
            restriction = [[sum(dot(a, v) << i for i, v in enumerate(basis))
                            for a in variables]]
            # The temporary nested form keeps the hot loop simple below.
            pattern = restriction[0]
            for ell in range(1, 1 << d):
                kval = 0
                for x in reps:
                    kval = max(kval, sum(pattern[labels[x ^ v] - 1] == ell
                                           for v in elements))
                coeff = [-1] * len(variables)
                for a, p in zip(variables, pattern):
                    if p == ell:
                        coeff[index[a]] += h
                rhs = 4 * kval - 2 * h
                rows.append(coeff)
                metadata.append((basis, ell, kval, h))
    return variables, rows, metadata


def solve_float(A, rhs, lower):
    import numpy as np
    from scipy.optimize import linprog
    # A W >= rhs, W >= lower; scipy takes -A W <= -rhs.
    result = linprog(np.ones(len(A[0])), A_ub=-np.asarray(A, float),
                     b_ub=-np.asarray(rhs, float),
                     bounds=[(lower, None)] * len(A[0]), method="highs")
    if not result.success:
        raise RuntimeError(result.message)
    return result


def gauss_solve(rows, rhs):
    """Return one exact solution of a full-column-rank square system."""
    n = len(rows[0])
    a = [[Fraction(x) for x in row] + [Fraction(y)]
         for row, y in zip(rows, rhs)]
    pivot = 0
    for col in range(n):
        q = next((r for r in range(pivot, len(a)) if a[r][col]), None)
        if q is None:
            continue
        a[pivot], a[q] = a[q], a[pivot]
        scale = a[pivot][col]
        a[pivot] = [x / scale for x in a[pivot]]
        for r in range(len(a)):
            if r != pivot and a[r][col]:
                scale = a[r][col]
                a[r] = [x - scale * y for x, y in zip(a[r], a[pivot])]
        pivot += 1
        if pivot == n:
            break
    if pivot < n:
        return None
    return [a[i][-1] for i in range(n)]


def gauss_rank(rows):
    a = [[Fraction(x) for x in row] for row in rows]
    rank = 0
    n = len(a[0]) if a else 0
    for col in range(n):
        q = next((r for r in range(rank, len(a)) if a[r][col]), None)
        if q is None:
            continue
        a[rank], a[q] = a[q], a[rank]
        scale = a[rank][col]
        a[rank] = [x / scale for x in a[rank]]
        for r in range(rank + 1, len(a)):
            if a[r][col]:
                scale = a[r][col]
                a[r] = [x - scale * z for x, z in zip(a[r], a[rank])]
        rank += 1
    return rank


def rational_primal(A, rhs, result, lower):
    n = len(A[0])
    x = result.x
    active = [(list(A[i]), rhs[i]) for i, v in enumerate(result.ineqlin.residual)
              if v < 2e-7]
    active += [([1 if j == i else 0 for j in range(n)], lower)
               for i, v in enumerate(result.lower.residual) if v < 2e-7]
    chosen = []
    rank = 0
    for row, y in active:
        newrank = gauss_rank([r for r, _ in chosen] + [row])
        if newrank > rank:
            chosen.append((row, y))
            rank = newrank
            if rank == n:
                break
    if rank < n:
        raise RuntimeError("active set did not span primal variables")
    qx = gauss_solve([r for r, _ in chosen], [y for _, y in chosen])
    return qx


def rational_dual(A, rhs, result):
    """Rationalize scipy's dual, then check y,z >= 0 and stationarity."""
    n = len(A[0])
    y = [Fraction(float(v)).limit_denominator(1000000)
         for v in (-result.ineqlin.marginals)]
    # For a minimization with lower bounds, z is the lower-bound dual.
    z = [Fraction(1) - sum(Fraction(A[i][j]) * y[i]
                           for i in range(len(A))) for j in range(n)]
    return y, z


def exact_check(A, rhs, primal, dual, lower):
    y, z = dual
    assert len(primal) == len(z) == len(A[0])
    assert len(y) == len(A) == len(rhs)
    assert all(x >= lower for x in primal)
    denominator = 1
    for x in primal:
        denominator = denominator // gcd(denominator, x.denominator) * x.denominator
    ints = [x.numerator * (denominator // x.denominator) for x in primal]
    slacks = [sum(c * x for c, x in zip(row, ints)) - b * denominator
              for row, b in zip(A, rhs)]
    assert min(slacks) >= 0, min(slacks)
    assert all(v >= 0 for v in y), min(y)
    assert all(v >= 0 for v in z), min(z)
    support = [i for i in range(len(y)) if y[i]]
    stationarity = [sum(A[i][j] * y[i] for i in support) + z[j]
                    for j in range(len(primal))]
    assert stationarity == [Fraction(1)] * len(primal), stationarity
    primal_obj = sum(primal)
    dual_obj = sum(Fraction(b) * yi for b, yi in zip(rhs, y)) + lower * sum(z)
    assert primal_obj == dual_obj, (primal_obj, dual_obj)
    # Complementarity is a useful guard against an accidentally unrelated
    # rational dual pair.
    assert all(si * yi == 0 for si, yi in zip(slacks, y))
    assert all((x - lower) * zi == 0 for x, zi in zip(primal, z))
    return (primal_obj, Fraction(min(slacks), denominator),
            len(support), sum(zi != 0 for zi in z))


def parse_fraction(value):
    return Fraction(value)


def verify_fixture(path=None):
    """Verify checked primal/dual certificates without NumPy or SciPy."""
    if path is None:
        path = Path(__file__).with_name("walsh_simultaneous_coset_weight_lp_20260915_fixture.json")
    with open(path, encoding="utf-8") as handle:
        fixture = json.load(handle)
    for name, spec in fixture.items():
        if name == "m32_archived":
            labels = list(M32)
        elif name == "m64_random_20260915":
            labels = map_random(64, 20260915)
        else:
            raise ValueError("unknown fixture " + name)
        variables, A, metadata = constraints(labels)
        rhs = [4 * m[2] - 2 * m[3] for m in metadata]
        primal = [parse_fraction(x) for x in spec["primal"]]
        if name == "m64_random_20260915":
            assert primal == [Fraction(63, 8) + Fraction((-1) ** dot(4, a)
                                                       + (-1) ** dot(29, a), 16)
                              for a in variables]
            counts = Counter(labels)
            q = [primal[a - 1] - Fraction(counts[a], 16) for a in variables]
            for d in range(1, 64):
                value = sum(x * (-1) ** dot(a, d) for a, x in zip(variables, q))
                expected = (Fraction(-4) if d == 4 else Fraction(-31, 8) if d == 29
                            else Fraction(-63, 8) if d % 2 else Fraction(-8))
                assert value == expected
        y = [Fraction(0)] * len(A)
        for i, value in spec.get("dual", {}).items():
            y[int(i)] = parse_fraction(value)
        z = [Fraction(0)] * len(variables)
        for i, value in spec.get("lower_dual", {}).items():
            z[int(i)] = parse_fraction(value)
        optimum, min_slack, ny, nz = exact_check(A, rhs, primal, (y, z), 5)
        assert optimum == parse_fraction(spec["objective"])
        print(name, "verified exact W=", optimum, "constraints=", len(A),
              "dual_support=", ny, nz, "min_slack=", min_slack)


def run(name, labels, lower=5):
    variables, A, metadata = constraints(labels)
    rhs = [4 * m[2] - 2 * m[3] for m in metadata]
    result = solve_float(A, rhs, lower)
    # The float solution is used only to find the active basis.  Fraction
    # reconstruction is deliberately bounded-denominator and exact-checked.
    primal = rational_primal(A, rhs, result, lower)
    dual = rational_dual(A, rhs, result)
    optimum, min_slack, ny, nz = exact_check(A, rhs, primal, dual, lower)
    gains = Counter(2 * m[2] - m[3] for m in metadata)
    print(name, "vars", len(variables), "constraints", len(A),
          "float", result.fun, "exact", optimum,
          "ratio", float(optimum) / len(variables),
          "max_gain", max(gains), "gain_hist", dict(sorted(gains.items())),
          "dual_support", ny, nz, "min_slack", min_slack)
    print("  primal", [str(x) for x in primal])
    return optimum, primal, dual


if __name__ == "__main__":
    if "--solve" in sys.argv:
        # PYTHONPATH may point at the unpacked cached SciPy wheel.
        run("M32 archived", list(M32))
        run("M64 random-20260915", map_random(64, 20260915))
    else:
        verify_fixture()
