"""Exact integer certificates for angles on multiplicatively modified Gaussian ideals.

Run with python3 graft_check.py. The three-ray search was numerical; this
certificate checks only the displayed candidate using exact arithmetic.
"""
from fractions import Fraction
from itertools import combinations
from math import gcd


def norm(v):
    return sum(t * t for t in v)


def det(v, w):
    return v[0] * w[1] - v[1] * w[0]


def mul(v, w):
    return v[0] * w[0] - v[1] * w[1], v[1] * w[0] + v[0] * w[1]


def qadd(v, w):
    return v[0] + w[0], v[1] + w[1]


def qscale(k, v):
    return k * v[0], k * v[1]


def qmul(v, w):
    return v[0] * w[0] + 2 * v[1] * w[1], v[0] * w[1] + v[1] * w[0]


def qpositive(v):
    """Whether a + b sqrt(2) is strictly positive."""
    a, b = v
    if a >= 0 and b >= 0:
        return a + b > 0
    if a <= 0 and b <= 0:
        return False
    return a * a > 2 * b * b if a > 0 else 2 * b * b > a * a


def three_fibres():
    beta = [(533, 8422), (1613, 994), (139, 326)]
    heights = [29**j * norm(vector) for j, vector in enumerate(beta)]
    assert heights == [71214173, 104104345, 105627077]
    assert all(gcd(*v) == 1 and norm(v) % 29 for v in beta)
    # Real and imaginary parts in Z[sqrt(2)] of beta_j (1+i sqrt(2))^j.
    # Multiply the third ray by the Gaussian unit -i to enter quadrant I.
    (a0,b0), (a1,b1), (a2,b2) = beta
    rays = [((a0, 0), (b0, 0)), ((a1, -b1), (b1, a1)),
            ((-b2, 2*a2), (a2, 2*b2))]
    assert all(qpositive(x) and qpositive(y) for x, y in rays)
    X = max(heights)
    for (x,y), (u,v) in combinations([rays[1],rays[2],rays[0]], 2):
        area = qadd(qmul(x,v), qscale(-1, qmul(y,u)))
        dot = qadd(qmul(x,u), qmul(y,v))
        assert qpositive(area)
        # 0 < angle <= tan(angle) = area / dot < 0.0216 / X.
        assert qpositive(qadd(qscale(216,dot), qscale(-10000*X,area)))
    print(f"p=29, tan(phi)=sqrt(2), heights {heights}, X*angular diameter < 0.0216")


def integer_bezout(a, b):
    """Return s,t with a*s+b*t=1."""
    r, rr, s, ss, t, tt = a, b, 1, 0, 0, 1
    while rr:
        q = r // rr
        r, rr, s, ss, t, tt = rr, r - q*rr, ss, s - q*ss, tt, t - q*tt
    assert abs(r) == 1
    return r*s, r*t


def rational_fibres(base, k, p=5):
    direction = (1, 0)
    for _ in range(k):
        direction = mul(direction,base)
    Q, h = norm(direction), p**k
    a, b = direction
    s, t = integer_bezout(a, b)
    alpha0 = (-t, s)
    while norm(alpha0) % p == 0:
        alpha0 = (alpha0[0]+a, alpha0[1]+b)
    assert det(direction, alpha0) == 1
    # The rational covector ell satisfies ell(alpha0)/ell(direction)=sqrt(h/Q).
    # With even k all coefficients can be taken integral.
    num, den = p**(k//2), norm(base)**(k//2)
    kernel = (den*alpha0[0]-num*a, den*alpha0[1]-num*b)
    ell = (-kernel[1], kernel[0])
    def parabolic(v, n):
        u = n*p*Q*(ell[0]*v[0]+ell[1]*v[1])
        return v[0]+u*kernel[0], v[1]+u*kernel[1]
    for n in range(1,101):
        alpha = parabolic(alpha0,n)
        Bbeta = parabolic(direction,n)
        numerator = mul((a,-b),Bbeta)
        assert all(x % Q == 0 for x in numerator)
        beta = tuple(x//Q for x in numerator)
        assert mul(direction,beta) == Bbeta
        assert det(Bbeta,alpha) == 1
        assert gcd(*alpha)==1 and gcd(*beta)==1 and (norm(alpha)*norm(beta))%p
        X = max(h*norm(beta),norm(alpha))
        normalized_sq = Fraction(X*X, norm(Bbeta)*norm(alpha))
        if normalized_sq <= Fraction(102,100)*Fraction(h,Q):
            print(f"base={base}, k={k}, n={n}, X has {len(str(X))} decimal digits,"
                  f" 1 <= normalized gap squared / (p/Norm(base))^k <= 1.02")
            assert normalized_sq >= Fraction(h,Q)
            return
    raise AssertionError("no witness")


if __name__ == "__main__":
    three_fibres()
    for base in [(3,2),(3,4)]:
        for k in (2,4,8):
            rational_fibres(base,k)
