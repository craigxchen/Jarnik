"""Finite checks for the deleted-row Hadamard short-character argument.

The endpoint and linear-form estimates in the companion note are theorem
inputs.  This checker audits the exact finite part: restriction of a
normalized Sylvester Hadamard matrix, whole cut classes, the Parseval good
old-column bound, constructive subset collisions, and the resulting row
certificate.  In particular, deleted-row old-column sums are included in
the collision key.  A retained-only collision is exhibited as a deliberate
failure case.
"""

from itertools import combinations
from math import comb, isqrt, log, sqrt
from random import Random


K = 2**40


def sylvester(order):
    """Return the normalized Sylvester Hadamard matrix of order 2^r."""
    assert order >= 1 and order & (order - 1) == 0
    rows = [[1]]
    while len(rows) < order:
        rows = [row + row for row in rows] + [row + [-x for x in row]
                                               for row in rows]
    assert len(rows) == order
    assert all(sum(rows[i][j] * rows[i][h] for i in range(order))
               == (order if j == h else 0)
               for j in range(order) for h in range(order))
    return rows


def column(matrix, index, rows=None):
    if rows is None:
        rows = range(len(matrix))
    return tuple(matrix[i][index] for i in rows)


def canonical(sign_column):
    """Canonical representative of an unoriented sign cut."""
    assert sign_column and all(value in (-1, 1) for value in sign_column)
    return tuple(value * sign_column[0] for value in sign_column)


def audit_cut_classes(matrix, deleted, extras):
    """Check nonconstant, distinct whole cut classes on both row sets."""
    n = len(matrix)
    retained = [i for i in range(n) if i not in set(deleted)]
    all_columns = [column(matrix, j) for j in range(1, n)] + list(extras)
    retained_columns = [column(matrix, j, retained) for j in range(1, n)]
    retained_columns += [tuple(extra[i] for i in retained) for extra in extras]
    for columns in (all_columns, retained_columns):
        assert all(len(set(values)) == 2 for values in columns)
        classes = [canonical(values) for values in columns]
        assert len(set(classes)) == len(classes)
    # The old core is already a valid restricted family when d<N/2.  Keep
    # this explicit check close to the theorem's hypothesis.
    assert len(deleted) * 2 < n
    old = retained_columns[: n - 1]
    assert len(set(map(canonical, old))) == n - 1


def make_extras(matrix, deleted, count, seed):
    """Produce arbitrary sign columns avoiding all old/extra cut classes."""
    rng = Random(seed)
    n = len(matrix)
    deleted = set(deleted)
    retained = [i for i in range(n) if i not in deleted]
    full_classes = {canonical(column(matrix, j)) for j in range(1, n)}
    retained_classes = {canonical(column(matrix, j, retained))
                        for j in range(1, n)}
    extras = []
    while len(extras) < count:
        candidate = tuple(1 if rng.getrandbits(1) else -1 for _ in range(n))
        restricted = tuple(candidate[i] for i in retained)
        if len(set(candidate)) < 2 or len(set(restricted)) < 2:
            continue
        if canonical(candidate) in full_classes:
            continue
        if canonical(restricted) in retained_classes:
            continue
        extras.append(candidate)
        full_classes.add(canonical(candidate))
        retained_classes.add(canonical(restricted))
    return extras


def parseval_good_columns(matrix, deleted, extras):
    """Return retained old-column Fourier sums and the Parseval good set."""
    assert extras
    n = len(matrix)
    retained = [i for i in range(n) if i not in set(deleted)]
    k = len(extras)
    numerators = [[sum(matrix[i][j + 1] * extra[i] for i in retained)
                   for extra in extras]
                  for j in range(n - 1)]
    total = sum(value * value for row in numerators for value in row)
    # Extend each extra by zero on deleted rows.  Parseval holds in the
    # full Hadamard basis; the restricted columns need not be orthogonal.
    for ell, extra in enumerate(extras):
        constant = sum(extra[i] for i in retained)
        assert constant * constant + sum(row[ell] ** 2 for row in numerators) == n * len(retained)
    assert total <= k * n * len(retained) <= k * n * n
    threshold_numerator = 2 * k * n * n
    good = [j for j, row in enumerate(numerators)
            if (n - 1) * sum(value * value for value in row)
            <= threshold_numerator]
    assert len(good) >= (n - 1) // 2
    bound = isqrt(threshold_numerator // (n - 1))
    assert all(abs(value) <= bound for j in good for value in numerators[j])
    return numerators, good, bound


def collision_key(matrix, deleted, numerators, subset, include_deleted=True):
    """Key of a q-subset, including deleted old-column sums when requested."""
    key = tuple(sum(numerators[j][ell] for j in subset)
               for ell in range(len(numerators[0]) if numerators else 0))
    if include_deleted:
        key += tuple(sum(matrix[row][j + 1] for j in subset)
                     for row in deleted)
    return key


def subset_collision(matrix, deleted, numerators, candidates, q,
                     include_deleted=True, max_subsets=100_000):
    """Construct the first equal-q-subset collision, if one is enumerable."""
    if q > len(candidates):
        return None
    count = comb(len(candidates), q)
    if count > max_subsets:
        return None
    seen = {}
    for subset in combinations(candidates, q):
        key = collision_key(matrix, deleted, numerators, subset,
                            include_deleted=include_deleted)
        if key in seen:
            prior = seen[key]
            v = [0] * (len(matrix) - 1)
            for j in subset:
                v[j] += 1
            for j in prior:
                v[j] -= 1
            assert any(v) and set(v) <= {-1, 0, 1}
            assert sum(abs(value) for value in v) <= 2 * q
            return v
        seen[key] = subset
    return None


def counting_criterion(order, deleted_count, extra_count, q):
    """The exact finite box inequality used by the pigeonhole step."""
    assert extra_count >= 1
    good_count = (order - 1) // 2
    bound = isqrt((2 * extra_count * order * order) // (order - 1))
    left = comb(good_count, q)
    right = (2 * q * bound + 1) ** extra_count * (2 * q + 1) ** deleted_count
    return left > right


def verify_certificate(matrix, deleted, extras, v):
    """Verify every coordinate and norm identity for the extracted character."""
    n = len(matrix)
    deleted = set(deleted)
    retained = [i for i in range(n) if i not in deleted]
    assert len(v) == n - 1 and any(v)
    assert set(v) <= {-1, 0, 1}
    ell = sum(abs(value) for value in v)
    lam_full = [sum(matrix[row][j + 1] * v[j] for j in range(n - 1))
                for row in range(n)]
    lam = [lam_full[row] for row in retained]

    # Constant-column orthogonality and all deleted coordinates are exact.
    assert sum(lam_full) == 0
    assert sum(lam) == 0
    assert all(lam_full[row] == 0 for row in deleted)
    assert sum(value * value for value in lam_full) == n * ell
    assert sum(value * value for value in lam) == n * ell
    l1 = sum(abs(value) for value in lam)
    assert l1 * l1 <= len(retained) * n * ell

    old_signature = [sum(lam[row_index] * matrix[row][j + 1]
                          for row_index, row in enumerate(retained))
                     for j in range(n - 1)]
    extra_signature = [sum(lam[row_index] * extra[row]
                           for row_index, row in enumerate(retained))
                       for extra in extras]
    assert old_signature == [n * value for value in v]
    assert extra_signature == [0] * len(extras)
    return ell, l1


def find_small_collision(matrix, deleted, extras, q_start, q_stop):
    """Search a few small q values; never enumerate a large subset family."""
    numerators, good, bound = parseval_good_columns(matrix, deleted, extras)
    del bound  # The caller separately audits the exact bound and criterion.
    for q in range(q_start, min(q_stop, len(good)) + 1):
        value = subset_collision(matrix, deleted, numerators, good, q)
        if value is not None:
            return value, q
    return None, None


def audit_fixture(order, deleted, extra_count, seed):
    matrix = sylvester(order)
    extras = make_extras(matrix, deleted, extra_count, seed)
    audit_cut_classes(matrix, deleted, extras)
    if not extras:
        return 0, 0

    numerators, good, bound = parseval_good_columns(matrix, deleted, extras)
    q = extra_count // 2 + 1
    criterion = counting_criterion(order, len(deleted), extra_count, q)
    if criterion:
        assert comb(len(good), q) > (
            2 * q * bound + 1) ** extra_count * (2 * q + 1) ** len(deleted)
    value, found_q = find_small_collision(matrix, deleted, extras, q,
                                           min(q + 2, 3))
    if value is not None:
        ell, l1 = verify_certificate(matrix, deleted, extras, value)
        assert ell <= 2 * found_q
        return 1, ell
    return 0, 0


def audit_k_zero(order, deleted):
    """For k=0, use a q=1 collision of deleted-row old-column patterns."""
    matrix = sylvester(order)
    extras = []
    audit_cut_classes(matrix, deleted, extras)
    numerators = [[] for _ in range(order - 1)]
    candidates = list(range(order - 1))
    value = subset_collision(matrix, deleted, numerators, candidates, 1)
    criterion = order - 1 > 2 ** len(deleted)
    if criterion:
        assert value is not None
    if value is None:
        return 0
    ell, l1 = verify_certificate(matrix, deleted, extras, value)
    assert ell == 2 and l1 * l1 <= (order - len(deleted)) * order * ell
    return 1


def deleted_constraint_regression():
    """A retained-only equal-sum collision can violate every deleted row."""
    order = 8
    deleted = (1, 2)
    matrix = sylvester(order)
    # This is nonconstant on both the full and retained row sets and is a
    # distinct cut class.  The retained-only q=1 collision is old columns
    # 1 and 2, while their deleted-row signs are opposite.
    extras = [(-1, -1, -1, -1, -1, -1, -1, 1)]
    audit_cut_classes(matrix, deleted, extras)
    numerators, good, _ = parseval_good_columns(matrix, deleted, extras)
    retained_only = subset_collision(matrix, deleted, numerators, good, 1,
                                     include_deleted=False)
    assert retained_only is not None
    assert all(sum(matrix[row][j + 1] * retained_only[j]
                   for j in range(order - 1)) != 0 for row in deleted)
    # Adding the d deleted-row coordinates forces q=2 here; the q=1
    # retained-only answer is therefore a deliberate false certificate.
    proper = subset_collision(matrix, deleted, numerators, good, 2,
                              include_deleted=True)
    assert proper is not None
    ell, _ = verify_certificate(matrix, deleted, extras, proper)
    assert ell <= 4
    return retained_only, proper


def counting_threshold_audit():
    """Check eventual fixed-(d,k) growth without enumerating subsets."""
    checked = 0
    for extra_count in range(1, 5):
        q = extra_count // 2 + 1
        for deleted_count in range(4):
            order = 8
            while (2 * deleted_count >= order or
                   not counting_criterion(order, deleted_count,
                                          extra_count, q)):
                order *= 2
                # The constants are deliberately crude: k=3,4 first cross
                # the exact box bound only at larger powers of two.  We
                # still evaluate an integer binomial, never its subsets.
                assert order < 2**45
            # A few later orders guard the claimed eventual direction.
            for _ in range(3):
                assert counting_criterion(order, deleted_count,
                                          extra_count, q)
                order *= 2
            checked += 1
    return checked


def k_zero_threshold_audit():
    checked = 0
    for deleted_count in range(7):
        order = 8
        while 2 * deleted_count >= order or order - 1 <= 2**deleted_count:
            order *= 2
        assert order - 1 > 2**deleted_count
        checked += 1
    return checked


def effective_arithmetic_audit():
    """Check Sol's explicit phase-gap inequalities at a safe large scale."""
    n = 10**18
    d = 3
    m = n - d
    for ell in (1, 2, 8, 32):
        logarithm = log(2 * n * sqrt(ell))
        assert m >= 16 * K * ell * (1 + logarithm)
        assert (m - 1) * log(5) > 8 * logarithm
    # The displayed 2N factor is the common-literal-unit C<=2 bound; it
    # also dominates the arbitrary-unit C<=sqrt(2) regime.


def explicit_cutoff(k, d):
    """Return the integer U(k,d) displayed in the uniform-bound note."""
    q = k // 2 + 1
    if k == 0:
        A = 2**d + 2
    else:
        A = ((4 * q)**(2 * q) * (25 * q * q * k)**k
             * (2 * q + 1)**(2 * d) + 1)
    return max(2 * d + 1, 4 * q, A, (256 * K * q)**2 + 1)


def explicit_cutoff_integer_audit():
    """Check (16)'s combinatorial cutoff using integer arithmetic only."""
    checked = 0
    for k in range(6):
        q = k // 2 + 1
        for d in range(5):
            U = explicit_cutoff(k, d)
            assert U >= 2 * d + 1 and U >= 4 * q
            if k == 0:
                # The q=1 deleted-pattern pigeonhole condition is strict.
                assert U - 1 > 2**d
            else:
                n = (U - 1) // 2
                B = isqrt((2 * k * U * U) // (U - 1))
                left = comb(n, q)
                right = ((2 * q * B + 1)**k
                         * (2 * q + 1)**d)
                assert left > right
                # Also verify the proof's integer lower-bound route directly.
                assert n >= U // 4
                assert B <= 2 * isqrt(k * U) + 1
                exponent = 2 * q - k
                constant = ((4 * q)**(2 * q) * (25 * q * q * k)**k
                            * (2 * q + 1)**(2 * d))
                assert U**exponent > constant
            checked += 1
    return checked


def explicit_cutoff_analytic_audit():
    """Numerically check the conservative inequalities following (16)."""
    checked = 0
    for k in range(6):
        q = k // 2 + 1
        ell = 2 * q
        for d in range(5):
            n = explicit_cutoff(k, d)
            m = n - d
            # Write log(e x) as 1+log(x), avoiding a rounded e constant.
            phase_log = log(2) + 1 + log(n) + 0.5 * log(ell)
            source_log = log(2) + log(n) + 0.5 * log(ell)
            assert m >= 16 * K * ell * phase_log
            assert (m - 1) * log(5) >= 8 * source_log
            # These are the elementary inequalities used to derive those
            # two estimates from N>(256Kq)^2 and N>1024.
            assert n > (256 * K * q)**2
            assert n > 1024
            assert 64 * K * q * log(n) <= n / 4
            assert 16 * log(n) <= n / 2
            checked += 1
    return checked


def family_audit():
    fixtures = ((8, ()), (8, (0,)), (8, (1, 2)), (8, (2, 5)),
                (16, (0, 3, 14)), (16, (1, 7, 12)),
                (32, (0, 5, 17, 30)), (64, ()))
    found = 0
    for order, deleted in fixtures:
        found += audit_k_zero(order, deleted)
        for extra_count in (1, 2, 3):
            hits, _ = audit_fixture(order, deleted, extra_count,
                                    seed=1000 + order + len(deleted) + extra_count)
            found += hits
    assert found >= 10
    return found


if __name__ == "__main__":
    false_v, proper_v = deleted_constraint_regression()
    found = family_audit()
    thresholds = counting_threshold_audit()
    zero_thresholds = k_zero_threshold_audit()
    effective_arithmetic_audit()
    integer_cutoffs = explicit_cutoff_integer_audit()
    analytic_cutoffs = explicit_cutoff_analytic_audit()
    print(f"PASS: {found} constructive Sylvester certificates across deleted-row fixtures.")
    print("PASS: exact cut classes, Parseval good-column bounds, all deleted coordinates, "
          "Hadamard signatures, and L1^2 <= MN ell.")
    print(f"PASS: retained-only false collision {false_v} rejected; proper collision {proper_v} verified.")
    print(f"PASS: {thresholds} fixed-(d,k) counting thresholds and {zero_thresholds} k=0 thresholds.")
    print(f"PASS: {integer_cutoffs} integer U(k,d) combinatorial cutoffs and "
          f"{analytic_cutoffs} conservative analytic cutoff checks.")
    print("PASS: Sol's K=2^40 arithmetic checked at N=10^18; no endpoint examples claimed.")
