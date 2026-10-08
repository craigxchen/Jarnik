"""Exact arithmetic checks of the median rotation normal form (stdlib)."""
from itertools import combinations, product
from math import gcd
from random import Random


def mul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def conj(a):
    return a[0], -a[1]


def norm(a):
    return a[0]**2+a[1]**2


def det(a, b):
    return a[0]*b[1]-a[1]*b[0]


def ggcd(a, b):
    while b != (0, 0):
        c = mul(a, conj(b))
        d = norm(b)
        q = ((2*c[0]+d)//(2*d), (2*c[1]+d)//(2*d))
        qb = mul(q, b)
        a, b = b, (a[0]-qb[0], a[1]-qb[1])
    return a


def primitive(a):
    d = gcd(*a)
    return a[0]//d, a[1]//d


def residue(a, b):
    n = norm(ggcd(a, b))
    assert det(a, b) % n == 0
    return abs(det(a, b))//n


def valuation(a, p):
    e = 0
    while a != (0, 0):
        z, n = mul(a, conj(p)), norm(p)
        if z[0] % n or z[1] % n:
            break
        a = z[0]//n, z[1]//n
        e += 1
    return e


def signed(a, p):
    return valuation(a, p)-valuation(a, conj(p))


def bezout(a, b):
    aa, bb, s, ss, t, tt = a, b, 1, 0, 0, 1
    while bb:
        q = aa//bb
        aa, bb = bb, aa-q*bb
        s, ss = ss, s-q*ss
        t, tt = tt, t-q*tt
    if aa == -1:
        aa, s, t = 1, -s, -t
    assert aa == 1 and s*a+t*b == 1
    return s, t


def maxnorm(a):
    return max(map(abs, a))


def check_rows(rows):
    p1, p2, p3 = rows[:3]
    g = ggcd(tuple(det(p2, p3)*x for x in p1),
             tuple(det(p1, p3)*x for x in p2))
    qs = [primitive(mul(p, conj(g))) for p in rows]
    assert all(norm(ggcd(q, conj(q))) == 1 for q in qs)
    assert all(norm(ggcd(a, b)) == 1 for a, b in combinations(qs[:3], 2))
    for i, j in combinations(range(len(rows)), 2):
        assert residue(rows[i], rows[j]) == residue(qs[i], qs[j])
    for pi in [(2, 1), (3, 2), (4, 1), (5, 2), (6, 1)]:
        assert signed(g, pi) == sorted(signed(p, pi) for p in rows[:3])[1]
        assert all(signed(q, pi) == signed(p, pi)-signed(g, pi)
                   for p, q in zip(rows, qs))
    x, y = qs[0]
    s, t = bezout(x, y)
    vs0 = [(s*u+t*v, -y*u+x*v) for u, v in qs]
    a0, b = vs0[1]
    # Search around the integer quotient to avoid rounding conventions.
    k = min(range(-a0//b-2, -a0//b+3), key=lambda j: abs(a0+j*b))
    vs = [(u+k*v, v) for u, v in vs0]
    a = vs[1][0]
    assert 2*abs(a) <= abs(b)
    T = max(residue(u, v) for u, v in combinations(qs, 2))
    assert vs[0] == (1, 0)
    assert max(maxnorm(v) for v in vs[:3]) <= 2*T
    col1 = qs[0]
    col2 = tuple((qs[1][j]-a*col1[j])//b for j in range(2))
    assert det(col1, col2) == 1
    for q, v in zip(qs, vs):
        assert q == tuple(col1[j]*v[0]+col2[j]*v[1] for j in range(2))
    H = max(maxnorm(q) for q in qs[:3])
    U = max(maxnorm(col1), maxnorm(col2))
    assert H <= 4*T*U and U <= 2*T*H


def main():
    rng = Random(91373)
    pool = [(x, y) for x in range(-25, 26) for y in range(-25, 26)
            if gcd(x, y) == 1 and (x+y) % 2]
    count = 0
    while count < 2500:
        rows = rng.sample(pool, 6)
        if any(det(a, b) == 0 for a, b in combinations(rows, 2)):
            continue
        check_rows(rows)
        count += 1
    print(f"Rotation, preserved residues, medians and small frame: {count} exact six-row cases pass.")
    for m in range(3, 13):
        cuts = list(product((0, 1), repeat=m))
        for i in range(m):
            count = sum(bits[i] != sorted(bits[:3])[1] for bits in cuts)
            assert count == 2**(m-2 if i < 3 else m-1)
        for i, j in combinations(range(m), 2):
            assert sum(b[i] != b[j] for b in cuts) == 2**(m-1)
    print("Full-profile median and separating-cut counts through 12 rows pass.")


if __name__ == "__main__":
    main()
