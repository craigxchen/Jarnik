"""Exact marks on q_P: A = prod_p f_p, f_p in {pi_p, conj pi_p}, unit eta = 1.
Own code; exact Python integers."""
import math, pickle, sys
from sympy import primerange

def two_sq(p):
    for a in range(1, math.isqrt(p) + 1):
        b2 = p - a * a
        b = math.isqrt(b2)
        if b * b == b2 and b > 0:
            return a, b
    raise ValueError

def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])

def build(P, cls=1, mod=4):
    ps = [p for p in primerange(P, 2 * P + 1) if p % mod == cls]
    pis = [two_sq(p) for p in ps]
    q = 1
    for p in ps:
        q *= p
    s = len(ps)
    marks = []
    for r in range(1 << s):
        z = (1, 0)
        for k, (a, b) in enumerate(pis):
            z = gmul(z, (a, b) if (r >> k) & 1 else (a, -b))
        assert z[0] ** 2 + z[1] ** 2 == q
        marks.append(z)
    return ps, pis, q, marks

if __name__ == "__main__":
    P = int(sys.argv[1])
    ps, pis, q, marks = build(P)
    with open(f"marks_{P}.pkl", "wb") as f:
        pickle.dump((ps, pis, q, marks), f)
    print(P, len(ps), q)
