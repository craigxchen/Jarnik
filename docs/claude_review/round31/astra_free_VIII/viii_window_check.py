"""Part VIII exact checks. Run python3 viii_window_check.py (requires SymPy)."""

from itertools import combinations
from math import lcm, prod

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form


def multiply(z, w):
    a, b = z
    c, d = w
    return a * c - b * d, a * d + b * c


def conjugate(z):
    return z[0], -z[1]


def norm(z):
    return sum(a * a for a in z)


def valuation(z, prime):
    p = norm(prime)
    count = 0
    while True:
        a, b = multiply(z, conjugate(prime))
        if a % p or b % p:
            return count
        z = a // p, b // p
        count += 1


def quadrant(z):
    x, y = z
    if x > 0 and y >= 0:
        return 0
    if x <= 0 and y > 0:
        return 1
    if x < 0 and y <= 0:
        return 2
    assert x >= 0 and y < 0
    return 3


def lifted_quarter_turns(primes, row):
    # Every multiplying factor has |arg| < pi/2. Its crossing of a quadrant
    # boundary is determined exactly by the quadrants before and after.
    z = 1, 0
    turns = 0
    for prime, power in zip(primes, row):
        factor = prime if power >= 0 else conjugate(prime)
        for _ in range(abs(power)):
            new = multiply(z, factor)
            step = (quadrant(new) - quadrant(z)) % 4
            assert step != 2
            turns += {0: 0, 1: 1, 3: -1}[step]
            z = new
    return turns, z


def check_marks():
    points = [(24791, 23712), (24745, 23760), (24636, 23873),
              (24567, 23944), (24260, 24255), (24255, 24260),
              (23944, 24567), (23873, 24636), (23760, 24745),
              (23712, 24791)]
    primes = [(2, 1), (3, 2), (4, 1), (5, 2), (8, 7)]
    exponents = [3, 2, 1, 1, 1]
    parent = prod(norm(prime) ** e for prime, e in zip(primes, exponents))
    assert parent == 1176852625 and all(norm(z) == parent for z in points)
    assert all(multiply(right, conjugate(left))[1] > 0
               for left, right in zip(points, points[1:]))
    x, y = points[0]
    assert points[-1] == (y, x)
    # Width is 2 atan((x-y)/(x+y)), and atan(u)<u for u>0.
    assert parent * (8 * (x - y)) ** 4 < (33 * (x + y)) ** 4
    A = [[valuation(z, prime) for prime in primes] for z in points]
    assert all(min(column) == 0 and max(column) == e
               for column, e in zip(zip(*A), exponents))
    B = Matrix([[2 * (v - base) for v, base in zip(row, A[0])]
                for row in A[1:]])
    k = []
    for row, point in zip(B.tolist(), points[1:]):
        turns, product = lifted_quarter_turns(primes, row)
        k.append(turns)
        remainder = multiply(product, [(1, 0), (0, -1), (-1, 0), (0, 1)][turns % 4])
        ratio = multiply(point, conjugate(points[0]))
        check = multiply(remainder, conjugate(ratio))
        assert check[0] > 0 and check[1] == 0
    smith = smith_normal_form(B, domain=ZZ)
    invariant_factors = [abs(smith[j, j]) for j in range(B.cols)]
    assert invariant_factors == [2, 2, 2, 2, 24]
    # H has columns forming a Z-basis of the column lattice of B-transpose.
    inverse_basis = hermite_normal_form(B.T).inv()
    orders = [lcm(*(int(x.q) for x in inverse_basis[:, j]))
              for j in range(B.cols)]
    assert orders == [24, 8, 6, 4, 24]
    exact_phase, parameters = B.gauss_jordan_solve(Matrix(k))
    assert parameters.rows == 0 and B * exact_phase == Matrix(k)
    assert k == [0, 3, 0, 1, 0, 1, -2, 1, 1]
    assert B.rank() == B.row_join(Matrix(k)).rank() == 5
    print("COMPUTED integer N =", parent, "; marks =", len(points), "; C < 33/4")
    print("ranks 5, 5; Smith invariants", invariant_factors, "; orders", orders)
    print("winding vector", k, "; exact phase in quarter turns", list(exact_phase))


    features = [row[:] for row in A]
    ranks = []
    for ell in (3, 7, 11, 19):
        colors = [(x % ell, y % ell) for x, y in points]
        used = sorted(set(colors))
        assert len(used) <= ell + 1
        for row, color in zip(features, colors):
            row.extend(int(color == c) for c in used)
        differences = Matrix([[v - w for v, w in zip(row, features[0])]
                              for row in features[1:]])
        ranks.append((ell, len(used), differences.rank()))
    assert ranks[-1][-1] == len(points) - 1
    print("(inert prime, observed classes, cumulative affine rank)", ranks)


def check_spread_background(q, K, S):
    # For P=q^(199/100), compare each distance to 10/sqrt(P) exactly.
    def near_center(numerator):
        remainder = numerator % (6*q)
        distance_units = min(remainder, 6*q - remainder)
        return distance_units**200 * q**199 <= (120*q)**200

    unsafe = [j for j in range(q) if any(
        near_center(q + 4*j + 4*q*value) for value in range(3))]
    selected = set(S)
    outside = [i for i in range(q) if i not in selected]
    u = [0] * q
    for i in S:
        total = sum(k*j for j, k in enumerate(K[i]))
        assert total % q == 0
        u[i] = total // q % 3

    # Solve (K^T u)_unsafe=0 by row reduction mod 3.
    rows = [[K[i][j] % 3 for i in outside] +
            [(-sum(K[i][j]*u[i] for i in S)) % 3] for j in unsafe]
    pivots = []
    for column in range(len(outside)):
        r = len(pivots)
        pivot = next((k for k in range(r, len(rows)) if rows[k][column]), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        factor = rows[r][column]  # nonzero elements of F_3 are self-inverse
        rows[r] = [factor*x % 3 for x in rows[r]]
        for k, row in enumerate(rows):
            multiplier = row[column]
            if k != r and multiplier:
                rows[k] = [(x - multiplier*y) % 3 for x, y in zip(row, rows[r])]
        pivots.append(column)
        if len(pivots) == len(rows):
            break
    assert len(pivots) == len(unsafe)
    for row, column in zip(rows, pivots):
        u[outside[column]] = row[-1]

    factors = [(-pow(q, -1, 3)*sum(K[i][j]*u[i] for i in range(q))) % 3
               for j in range(q)]
    assert all(factors[j] == 0 for j in unsafe)
    numerators = [q + 4*j + 4*q*factors[j] for j in range(q)]
    assert all(not near_center(n) for n in numerators)
    turns = [sum(K[i][j]*numerators[j] for j in range(1, q)) for i in S]
    assert all((value - turns[0]) % (6*q) == 0 for value in turns)
    for a, b in combinations(numerators[1:], 2):
        assert (a-b) % (12*q) not in (0, 1, 2, 3, 12*q-1, 12*q-2, 12*q-3)
        assert (a+b) % (12*q) not in (0, 1, 12*q-1)
    print("COMPUTED spread q =", q, "; prescribed labels =", len(unsafe),
          "; rank over F_3 =", len(pivots), "; integer row winding checks =", len(S))


def paley_moment_lifts(q=383):
    # Check the q/4 exact congruence moments after deleting column zero.
    # Interpolation gives a deterministic completion on the complement of S.
    squares = {i*i % q for i in range(1, q)}
    chi = [0 if i == 0 else (1 if i in squares else -1) for i in range(q)]
    K = [[chi[(j-i) % q] for j in range(q)] for i in range(q)]
    positive = [i for i in range(q) if K[i][0] == 1]
    assert len(positive) == (q-1)//2
    S = positive[:-2]
    check_spread_background(q, K, S)
    keep = range(1, q)
    distances = [sum(abs(K[i][j] - K[k][j]) for j in keep)
                 for i, k in combinations(S, 2)]
    assert set(distances) == {q+1}
    # Off the marked nodes evaluate each basis polynomial l_i as
    # product_k(j-k) divided by (j-i) times product_{k != i}(i-k).
    node_products = [prod(j-k for k in S) % q for j in range(q)]
    inverses = [0] + [pow(k, -1, q) for k in range(1, q)]
    interpolation = []
    for i in S:
        den_inv = inverses[prod(i-k for k in S if k != i) % q]
        interpolation.append([1 if i == j else den_inv * node_products[j] * inverses[(j-i) % q] % q
                              for j in range(q)])
    def centered(x):
        return (x + q//2) % q - q//2
    max_size = 0
    J = q//4
    for power in range(1, J+1):
        digit = [centered(pow(j, power, q)) for j in range(q)]
        value = [sum(a*b for a, b in zip(row, digit)) for row in K]
        assert all(x % q == 0 for x in value)
        tau = [x//q for x in value]
        completion = [centered(sum(tau[i]*basis[j] for i, basis in zip(S, interpolation)))
                      for j in range(q)]
        for i in S:
            completion[i] = tau[i]
        # K is skew, so d - K^T completion = d + K completion.
        lift = [d + sum(a*b for a, b in zip(row, completion))
                for row, d in zip(K, digit)]
        assert all((d-a) % q == 0 for d, a in zip(lift, digit))
        assert len({sum(K[i][j]*lift[j] for j in keep) for i in S}) == 1
        max_size = max(max_size, max(map(abs, lift)))
    print("COMPUTED Paley q =", q, "; marks =", len(S), "; columns =", len(keep),
          "; pairs =", len(distances), "; all pair lengths =", q+1,
          "; unweighted margin =", q+1-len(keep), "; integer moment lifts =", J,
          "; max lift sup-norm =", max_size)


if __name__ == "__main__":
    check_marks()
    paley_moment_lifts()
