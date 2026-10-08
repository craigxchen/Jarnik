"""Exact bounded checks for the prime-field polynomial character gap."""

from itertools import product
from math import comb, isqrt
from random import Random

from check_general_four_row_quartic_slack_obstruction import paley


def verify(matrix, chi, c, polynomial=False):
    m, q = len(c), len(c) - 1
    h = sum(x != 0 for x in c)
    assert h and sum(c) == 0
    transform = [sum(c[i] * matrix[i][j] for i in range(m)) for j in range(m)]
    l1 = sum(map(abs, transform))
    active = sum(t != 0 for t in transform)
    assert transform[0] == 0 and all(t % 2 == 0 for t in transform)
    assert sum(t*t for t in transform) == m*h
    assert l1 >= m
    assert 2*h*active <= (h+2)*l1 - m*h
    if h < q:
        pvalues = [sum(chi(x-y)*c[y+1] for y in range(q)) + c[0]
                   for x in range(q)]
        assert pvalues == [transform[x+1]+c[x+1] for x in range(q)]
        assert max(map(abs, pvalues)) <= h
        assert any(x % q for x in pvalues)
        required = m//2 + (c[0] == 0)
        assert sum(x % q != 0 for x in pvalues) >= required
        assert active+h >= m//2+1
        assert (h+2)*l1 >= 2*h*(m+1-h)
        if polynomial:
            degree = (q-1)//2
            coeff = [sum(c[y+1]*comb(degree,k)*pow(-y,degree-k,q)
                         for y in range(q)) % q for k in range(degree+1)]
            coeff[0] = (coeff[0]+c[0]) % q
            assert any(coeff)
            if c[0] == 0:
                assert coeff[-1] == 0
            assert all(sum(a*pow(x,k,q) for k,a in enumerate(coeff)) % q
                       == pvalues[x] % q for x in range(q))
    if m > 12 and l1 == m:
        assert h in (2, m//2, m)
    return transform, l1, active


def exhaustive(q):
    matrix, chi = paley(q)
    m, count = q+1, 0
    norm_one = set()
    for c in product((-1,0,1), repeat=m):
        if not any(c) or sum(c):
            continue
        _, l1, _ = verify(matrix, chi, c)
        if l1 == m:
            norm_one.add(c)
        count += 1
    expected = set()
    for i in range(m):
        for j in range(m):
            if i != j:
                expected.add(tuple(int(k==i)-int(k==j) for k in range(m)))
    for a in range(1,m):
        for sign in (-1,1):
            expected.add(tuple(sign*matrix[i][a] for i in range(m)))
        for b in range(a+1,m):
            for sa,sb in product((-1,1), repeat=2):
                expected.add(tuple((sa*matrix[i][a]+sb*matrix[i][b])//2
                                   for i in range(m)))
    assert norm_one == expected
    print(f"q={q}: {count} ternary characters; {len(norm_one)} norm-one classified")


def larger_fixtures():
    rng = Random(19741)
    for q in (19, 43, 59, 83):
        matrix, chi = paley(q)
        m, b = q+1, 5
        # Each original label occurs once, except a random duplicate.
        labels = list(range(1,m)) + [rng.randrange(1,m)]
        rng.shuffle(labels)
        # Sample distinct physical slots from two copies of each label.
        # Require an omitted label to cover nonsurjective assignments too.
        slots = [(label,copy) for label in range(1,m) for copy in range(2)]
        for _ in range(100):
            chosen_slots = rng.sample(slots,m)
            nonsurjective_labels = [label for label,_ in chosen_slots]
            if len(set(nonsurjective_labels)) < m-1:
                break
        else:
            raise AssertionError("failed to sample an assignment with an omitted label")
        assert len(set(chosen_slots)) == m
        assert max(nonsurjective_labels.count(label) for label in set(nonsurjective_labels)) <= 2
        assignments = (labels,nonsurjective_labels)
        count = 0
        for h in range(2,m+1,2):
            for trial in range(12):
                rows = rng.sample(range(m),h)
                c = [0]*m
                for k,i in enumerate(rows):
                    c[i] = 1 if k<h//2 else -1
                transform, l1, _ = verify(matrix,chi,c,polynomial=trial==0)
                for assignment in assignments:
                    physical = b*l1//2
                    for i,label in enumerate(assignment):
                        baseline = transform[label]//2
                        physical += abs(baseline-c[i]*matrix[i][label])-abs(baseline)
                    if m >= 36:
                        assert 2*physical > b*(m-1)
                count += 1
        print(f"q={q}: {count} exact polynomial fixtures, {2*count} assignment checks")
    # Integer checks supplement the analytic concavity proof, without sqrt rounding.
    for m in range(36,2001,4):
        for h in range(4,isqrt(4*m)+1,2):
            twice_f = 5*(h-2)*m-12*h*h+11*h+10
            assert twice_f > 0
    # The endpoint cubic at x=6 and its derivative are positive; derivative increases.
    assert 5*6**3-29*6**2+11*6+5 == 107
    assert 15*6**2-58*6+11 > 0 and 30*6-58 > 0


if __name__ == "__main__":
    exhaustive(7)
    exhaustive(11)
    larger_fixtures()
    print("PASS: prime-field support/norm gap and capacity-two half-height margin")
