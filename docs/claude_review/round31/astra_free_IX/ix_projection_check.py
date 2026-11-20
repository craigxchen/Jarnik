"""Exact checks for Part IX's projection bound. Run python3 ix_projection_check.py."""
from itertools import combinations
from math import isqrt, prod

MARKS = [(24791, 23712), (24745, 23760), (24636, 23873),
         (24567, 23944), (24260, 24255), (24255, 24260),
         (23944, 24567), (23873, 24636), (23760, 24745),
         (23712, 24791)]
FACTORS = [(2, 1), (3, 2), (4, 1), (5, 2), (8, 7)]
EXPONENTS = (3, 2, 1, 1, 1)


def multiply(z, w):
    x, y = z
    a, b = w
    return x*a - y*b, x*b + y*a


def conjugate(z):
    return z[0], -z[1]


def norm(z):
    return sum(value * value for value in z)


def valuation(z, prime):
    p = norm(prime)
    count = 0
    while True:
        a, b = multiply(z, conjugate(prime))
        if a % p or b % p:
            return count
        z = a//p, b//p
        count += 1


def primes_through(n):
    return [p for p in range(2,n+1) if all(p % k for k in range(2,isqrt(p)+1))]


def local_projection_check():
    # Count every possible projection on every residue circle, including bad
    # reduction and both local splitting types, with no height estimate.
    counts = [0,0]
    for p in primes_through(251)[1:]:
        squares = {x*x % p for x in range(p)}
        for d in range(p):
            count = sum((d-x*x) % p in squares for x in range(p))
            if d:
                expected = (p - (1 if p % 4 == 1 else -1)
                            + (2 if d in squares else 0))//2
            else:
                expected = p if p % 4 == 1 else 1
            assert count == expected
            counts[p % 4 == 3] += 1
    print("projection images at all residues and odd primes <= 251:",
          "split",counts[0],"inert",counts[1])


def record_check():
    parent = prod(norm(p)**e for p,e in zip(FACTORS,EXPONENTS))
    assert parent == 1176852625 and all(norm(z)==parent for z in MARKS)
    rows = [[valuation(z,p) for p in FACTORS] for z in MARKS]
    minima = [min(col) for col in zip(*rows)]
    maxima = [max(col) for col in zip(*rows)]
    assert minima == [0]*len(FACTORS) and tuple(maxima)==EXPONENTS
    certificates = []
    unit_factors = ((1,0),(0,1),(-1,0),(0,-1))
    for i,j in combinations(range(len(MARKS)),2):
        v = [x+y-e for x,y,e in zip(rows[i],rows[j],EXPONENTS)]
        g = (1,0)
        for p,t in zip(FACTORS,v):
            for _ in range(abs(t)):
                g = multiply(g, p if t>0 else conjugate(p))
        T = prod(norm(p)**abs(t) for p,t in zip(FACTORS,v))
        assert norm(g)==T
        unit = next(u for u in unit_factors if multiply(
            multiply(MARKS[i],MARKS[j]),conjugate(g)) ==
            multiply(multiply((parent,0),g),u))
        # Correct b/bar b by a unit, with at most one ramified factor.
        b = next(multiply(g,c) for c in ((1,0),(0,1),(1,1),(1,-1))
                 if c == multiply(unit,conjugate(c)))
        assert multiply(multiply(MARKS[i],MARKS[j]),conjugate(b)) == multiply((parent,0),b)
        assert norm(b) in (T,2*T)
        coords = [multiply(z,conjugate(b)) for z in MARKS]
        assert all(norm(w) == parent*norm(b) for w in coords)
        projections = sorted({x for x,y in coords})
        assert len(projections)>=len(MARKS)//2
        if T==1:
            certificates.append((i,j,b,projections))
    assert len(certificates)==5
    i,j,b,projections = certificates[0]
    assert b==(1,1) and projections==[48503,48505,48509,48511,48515]
    candidates = [x for x in range(projections[0],projections[-1]+1)
                  if isqrt(2*parent-x*x)**2 == 2*parent-x*x]
    assert candidates==projections
    print("ten-mark norm", parent, "pair certificates", 45, "T=1:", len(certificates),
          "b =",b,"projections",projections)


def ternary_witness():
    # The previous D=12 monic squarefree example. The imported routines use
    # python-flint only to represent F_9[t] and decode its base-nine encoding.
    from vii_squarefree_witness import decode, ring
    points = list(map(decode, [391358742215, 391358571511, 391360438826,
                               391357359206, 391358930855]))
    def bar(z):
        return ring([c**3 for c in z.coeffs()])
    parent = points[0]*bar(points[0])
    assert all(z*bar(z) == parent for z in points)
    histogram = {}
    for x,y in combinations(points, 2):
        product = x*y
        b = product // product.gcd(parent)
        assert product*bar(b) == parent*b
        traces = [z*bar(b) + bar(z)*b for z in points]
        span_degree = max((t-u).degree() for t,u in combinations(traces,2))
        assert span_degree == b.degree()
        assert len(set(map(str,traces))) == 4
        histogram[span_degree] = histogram.get(span_degree,0) + 1
    print("ternary D=12 five marks; pair normal degree equals trace diameter degree:",
          sorted(histogram.items()), "; four distinct trace values per pair")


if __name__ == '__main__':
    local_projection_check()
    record_check()
    ternary_witness()
