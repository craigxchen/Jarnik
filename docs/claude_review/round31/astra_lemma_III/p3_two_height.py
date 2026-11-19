"""Exact two-height exclusion for linear-block five-stars (Python and python-flint).

Polynomials are represented by sparse FLINT multivariates. Every discarded
irreducible factor divides a proved noncollision condition. Resultants are used
only as necessary conditions, including at parameter values where degrees drop.
"""
from itertools import combinations

from flint import fmpz_mpoly_ctx, fmpz_mod_mpoly_ctx


def determinant(rows):
    if len(rows) == 1:
        return rows[0][0]
    return sum((-1)**i * entry *
               determinant([r[:i] + r[i + 1:] for r in rows[1:]])
               for i, entry in enumerate(rows[0]))


def constraints(ctx):
    one = ctx.constant(1)
    t = [0*one, one, *ctx.gens()]
    z = [t[i] - t[(i + 1) % 5] + t[(i + 2) % 5] for i in range(5)]
    d = [(t[i] - t[(i + 1) % 5]) * (t[(i + 2) % 5] - t[(i + 1) % 5])
         for i in range(5)]
    f = [d[i] * t[(i + 1) % 5] for i in range(5)]

    def interpolate(values, count):
        matrix = [[z[i]**j for j in range(count)] for i in range(count)]
        denominator = determinant(matrix)
        coefficients = [determinant([row[:j] + [values[i]] + row[j + 1:]
                                     for i, row in enumerate(matrix)])
                        for j in range(count)]
        return coefficients, denominator

    q, qden = interpolate(d, 3)
    c, cden = interpolate(f, 4)
    equations = [sum(q[j]*z[i]**j for j in range(3)) - qden*d[i] for i in (3, 4)]
    equations += [sum(c[j]*z[4]**j for j in range(4)) - cden*f[4]]
    equations += [3*c[3]*qden**2 - cden*q[2]*(qden - 4*q[2]),
                  c[2]*qden**2 + 2*cden*q[1]*q[2]]
    excluded = [v[i] - v[j] for i, j in combinations(range(5), 2) for v in (t, z)]
    return [remove_collision_factors(p, excluded) for p in equations], excluded


def remove_collision_factors(p, excluded):
    assert p
    result = p.context().constant(1)
    for factor, _ in p.factor()[1]:
        if any(e.total_degree() == 1 and not e % factor for e in excluded):
            continue
        result *= factor
    return result


def check(ctx):
    x, _, _ = ctx.gens()
    equations, excluded = constraints(ctx)
    projections = [remove_collision_factors(equations[0].resultant(q, 2), excluded)
                   for q in equations[1:]]
    search = ctx.constant(1)
    for factor, _ in projections[0].factor()[1]:
        candidate = ctx.constant(0)
        for g in projections[1:]:
            r = factor.resultant(g, 1)
            if r:
                candidate = candidate.gcd(r) if candidate else r
        assert candidate
        search *= remove_collision_factors(candidate, excluded)
    p = remove_collision_factors(search, excluded)
    deg = p.degrees()[0]
    reverse = sum(int(c)*x**ex[0]*(x - 1)**(deg - ex[0])
                  for ex, c in p.to_dict().items())
    gcd = p.gcd(reverse)
    assert gcd.total_degree() == 1 and not gcd % (x - 2)
    print("eliminant degree", deg, "reversal gcd", gcd)


if __name__ == "__main__":
    names = ["x", "y", "z"]
    check(fmpz_mpoly_ctx.get(names))
    for characteristic in (17, 23, 103, 107):
        print("characteristic", characteristic)
        check(fmpz_mod_mpoly_ctx.get(names, characteristic))
