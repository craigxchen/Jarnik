"""Exact collapse and identity-agreement fixtures for growing-degree maps."""

from math import gcd


def multiply(p, q):
    out = [0]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return out


def add(p, q):
    assert len(p) == len(q)
    return [a+b for a, b in zip(p, q)]


def evaluate(p, x, t):
    d = len(p)-1
    return sum(c*x**(d-j)*t**j for j, c in enumerate(p))


def mul(z, w):
    return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])


def conj(z):
    return (z[0], -z[1])


def norm(z):
    return z[0]*z[0]+z[1]*z[1]


def nearest(n, d):
    return (n+d//2)//d if n >= 0 else -nearest(-n, d)


def divrem(z, w):
    n = mul(z, conj(w))
    q = (nearest(n[0], norm(w)), nearest(n[1], norm(w)))
    product = mul(q, w)
    return q, (z[0]-product[0], z[1]-product[1])


def gaussian_gcd(z, w):
    while w != (0, 0):
        _, remainder = divrem(z, w)
        z, w = w, remainder
    return z


def gaussian_lcm(z, w):
    g = gaussian_gcd(z, w)
    q, remainder = divrem(w, g)
    assert remainder == (0, 0)
    return mul(z, q)


def least_squared_radius(denominators):
    L = (1, 0)
    for denominator in denominators:
        L = gaussian_lcm(L, denominator)
    return norm(L)


def check_rows(rows, anchor_order):
    r = len(rows)
    factors = [[t, -x] for x, t in rows]
    product = [1]
    for factor in factors:
        product = multiply(product, factor)
    product_x = 1
    for x, _ in rows:
        product_x *= x
    assert product[-1] == (-1)**r*product_x

    # Collapse map P=X^D, Q=T^e prod_i(t_i X-x_i T).
    e = anchor_order
    D = r+e
    P = [1]+[0]*D
    Q = multiply([0]*e+[1], product)
    assert len(Q) == D+1 and Q[:e] == [0]*e
    H = max(max(map(abs, P)), max(map(abs, Q)))
    assert H >= product_x
    rho_abs = abs(Q[-1])**D  # Exact binary resultant Res(X^D,Q).
    assert rho_abs == product_x**D
    for x, t in rows:
        assert gcd(x, t) == 1 and t != 0
        p, q = evaluate(P, x, t), evaluate(Q, x, t)
        assert p == x**D and q == 0
        assert gcd(p, q) == x**D and rho_abs % gcd(p, q) == 0
    assert least_squared_radius([(1, 0)]*r) == 1

    # Nonidentity map agreeing with u=T/X at every selected row.
    D = r+e
    P = [1]+[0]*D
    identity_numerator = [0, 1]+[0]*(D-1)  # T X^(D-1)
    Q = add(identity_numerator, multiply([0]*e+[1], product))
    A = add(multiply([1, 0], Q),
            [-c for c in multiply([0, 1], P)])
    assert A == multiply([0]*e+[1, 0], product)  # XT^e prod factors.
    H = max(max(map(abs, P)), max(map(abs, Q)))
    assert max(map(abs, A)) <= 2*H
    assert max(map(abs, A)) >= product_x
    rho_abs = abs(Q[-1])**D
    source_denominators = []
    target_denominators = []
    for x, t in rows:
        p, q = evaluate(P, x, t), evaluate(Q, x, t)
        assert p == x**D and q == t*x**(D-1)
        content = gcd(p, q)
        assert content == x**(D-1) and rho_abs % content == 0
        primitive = (p//content, q//content)
        assert primitive == (x, t)
        source_denominators.append((x, -t))
        target_denominators.append(conj(primitive))
    assert least_squared_radius(source_denominators) == \
        least_squared_radius(target_denominators)
    assert least_squared_radius(source_denominators) > 1


def main():
    primes = [101, 103, 107, 109, 113, 127, 131]
    for r in range(2, 8):
        rows = [(primes[i], 2*(i+1)) for i in range(r)]
        for e in (1, 2, 3):
            check_rows(rows, e)
    print("PASS: 18 growing-degree exact collapse fixtures; coefficient and resultant bounds")
    print("PASS: 18 nonidentity maps fixing all selected rational directions and their exact least radii")


if __name__ == "__main__":
    main()
