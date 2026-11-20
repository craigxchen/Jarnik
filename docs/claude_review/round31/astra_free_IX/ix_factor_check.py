"""Exact certificates for the factored local witnesses in Part IX.
Run python3 ix_factor_check.py. Requires python-flint and SymPy.
"""
from itertools import combinations
from math import gcd, isqrt
from flint import fmpz_mat
from sympy import primerange

ONE = (1, 0)
I = (0, 1)


def conj(z):
    return z[0], -z[1]


def norm(z):
    return z[0]**2 + z[1]**2


def mul(x, y, b):
    return ((x[0]*y[0] - x[1]*y[1]) % b, (x[0]*y[1] + x[1]*y[0]) % b)


def power(z, exponent, b):
    if exponent < 0:
        inv = pow(norm(z), -1, b)
        z = mul((inv, 0), conj(z), b)
        exponent = -exponent
    result = ONE
    while exponent:
        if exponent % 2:
            result = mul(result, z, b)
        z = mul(z, z, b)
        exponent //= 2
    return result


def combined(values, exponents, b):
    ans = ONE
    for z, e in zip(values, exponents):
        if e:
            ans = mul(ans, power(z, int(e), b), b)
    return ans


def local_test():
    q = 11
    labels = list(range(1, q))
    sq = {i*i % q for i in labels}
    def K(i, j):
        return 0 if i == j else (1 if (j-i) % q in sq else -1)
    marks = [i for i in labels if K(i, 0) == 1]
    n, m = len(labels), len(marks)
    A = [[1+K(i, j) for j in labels] for i in marks]
    assert all({0, 2} <= set(column) for column in zip(*A))
    D = fmpz_mat([[y-x for x, y in zip(A[0], row)] for row in A[1:]])
    H = fmpz_mat([[int(i==j)-int(marks[0]==j) for j in labels] for i in marks[1:]])
    differences = H-D
    assert all(x % 2 == 0 for x in differences.entries())
    hermite, trans = D.transpose().hnf(transform=True)
    U = fmpz_mat([[trans[i, j] for i in range(m-1)] for j in range(n)])
    assert (D*U).is_one()
    V = U * fmpz_mat(m-1, n, [x//2 for x in differences.entries()])
    T = fmpz_mat([[int(i==j)+2*V[i,j] for j in range(n)] for i in range(n)])
    assert D*T == H
    # Deliberately include nonzero even row windings.
    winding = [2*i*i+4 for i in range(m)]
    unit_power = U * fmpz_mat([[(k-winding[0])//2] for k in winding[1:]])
    primes = [p for p in primerange(101, 500) if p % 4 == 1][:n]
    factors = []
    for p in primes:
        x = next(x for x in range(1, isqrt(p)+1) if isqrt(p-x*x)**2 == p-x*x)
        factors.append((x, isqrt(p-x*x)))

    def make_local(b, skip):
        gammas = [mul(z, power(conj(z), -1, b), b) if j != skip else ONE
                  for j, z in enumerate(factors)]
        local_marks = [power(I, -k, b) for k in winding]
        rr = []
        for j, z in enumerate(factors):
            if j == skip:
                rr.append(None)
                continue
            f = mul(z, power(I, int(unit_power[j,0]), b), b)
            f = mul(f, combined(gammas, [V[j,k] for k in range(n)], b), b)
            assert norm(f) % b == primes[j] % b
            rr.append(f)
            for u in range(m):
                e = A[u][j]
                local_marks[u] = mul(local_marks[u],
                    mul(power(f,e,b), power(conj(f),2-e,b), b), b)
        return gammas, local_marks, rr

    good_pairs = same_level = 0
    for ell in (3,5,7,11,13,17,19,23,29,31,37):
        for a in (1,2,3):
            b = ell**a
            chars, points, _ = make_local(b, None)
            for u,v in combinations(range(m), 2):
                assert mul(points[u], chars[marks[v]-1], b) == mul(points[v], chars[marks[u]-1], b)
                good_pairs += 1
    for j, p in enumerate(primes):
        for b in (p, p*p, p**3):
            chars, points, _ = make_local(b,j)
            for u,v in combinations(range(m), 2):
                if A[u][j] == A[v][j]:
                    assert j+1 not in (marks[u],marks[v])
                    assert mul(points[u], chars[marks[v]-1], b) == mul(points[v], chars[marks[u]-1], b)
                    same_level += 1
    print("q =",q,"marks =",marks,"maximum absolute right-inverse entry =",max(abs(x) for x in U.entries()))
    print("exterior pair identities =",good_pairs,"internal level identities =",same_level)
    # At a single external prime, violate the necessary height condition
    # for some singleton strengthened by every monomial's local collisions.
    for ell in primerange(19, 2000):
        if ell in primes:
            continue
        _, _, rr = make_local(ell,None)
        for j, r in enumerate(rr):
            ratio = mul(r, power(conj(r),-1,ell),ell)
            for e in range(4):
                if ratio == power(I,e,ell) and ell**2 > 2*primes[j]:
                    print("negative control: physical prime", primes[j],
                          "at external prime",ell,"ratio = i^",e,
                          "; squared divisor exceeds twice the weight prime")
                    return
    raise AssertionError("expected a failed singleton")


def real_power_test():
    # One exact coordinate-ray test for IX.7, with arbitrarily excluded
    # small primes modeled by this finite list.
    from math import lcm
    z = (2, 1)
    excluded = (2, 3, 5, 7, 13)
    period = 4
    for ell in excluded:
        if ell == norm(z):
            continue
        order = 1
        w = (z[0] % ell, z[1] % ell)
        while w != ONE:
            w = mul(w, z, ell)
            order += 1
        period = lcm(period, order)
    for k in range(1, 21):
        for ell in excluded:
            assert power(z, period*k, ell)[0] % ell
    # beta = arg(2+i). At assigned angle beta/3 and exponent 6 the
    # cosine has magnitude 3/5, versus the required 117/125. Its two
    # maximal exterior collision powers are 9 and 13 when E = {2,5}.
    assert mul(z, z, 1000) == (3,4)
    assert power(z, 6, 1000) == (883,44)  # real coordinate -117
    for modulus in (9, 13):
        ratio = mul(z, power(conj(z), -1, modulus), modulus)
        assert power(ratio, 6, modulus) == ((-1) % modulus, 0)
    assert 3*125 < 117*5
    print("clean real parts: period", period, "at twenty multiples, excluded", excluded,
          "; cube-root angle fails at exponent 6 and cutoff 13")


if __name__ == '__main__':
    local_test()
    real_power_test()
