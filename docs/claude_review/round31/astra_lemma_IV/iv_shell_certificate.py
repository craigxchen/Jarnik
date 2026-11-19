"""Exact polynomial certificates for both orientations in Part IV."""
from itertools import combinations
from math import prod, isqrt, gcd
from sympy import Poly, Rational, expand, isprime, symbols

P, q = symbols("P q")
edges = list(combinations(range(5), 2))
cycle = {(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)}

def check_balanced():
    shifts = (145, 79, -31, 109, 85, -47, 113, 19, -41, 25)
    J = {e: P + q*b for e, b in zip(edges, shifts)}
    assert len(set(shifts) | {0, 329}) == 12
    R = (4680, 3861, 2860, 3168, 2808)
    d = (-76, -55, -6, -34, -4)
    f = (26095, Rational(43223, 2), 16715, 16813, 16303)
    T = [P**2 + a*q*P + b*q**2 for a, b in zip(d, f)]
    incident = [{e for e in edges if (i in e) != (e in cycle)} for i in range(5)]
    n = [prod(J[e] for e in s) for s in incident]
    for i, j in edges:
        shared = incident[i] & incident[j]
        assert len(shared) == 2
        assert expand(T[i]*R[j]-T[j]*R[i]-(R[j]-R[i])*prod(J[e] for e in shared)) == 0
    for i in range(5):
        err = Poly(expand((P+329*q)*T[i]**2 - n[i]), P, q)
        assert all(a<=2 and a+b==5 for (a,b),_ in err.terms())
    print("Balanced: all ten signed determinants; all five degree-five shell identities.")

    # A squarefree specialization with private split blocks and p_* = K²+1.
    k, step = 3999318414, 924
    origin = k*k + 1 - 329*step
    factors = [
        [15994547776559305381], [15994547776559244397],
        [1669, 3595601, 2665287953], [9510509, 1681776209513],
        [41, 390110921379493901], [15994547776559127973],
        [15994547776559275813], [709, 22559305749730873],
        [15994547776559133517], [10513, 1521406618145077],
    ]
    primes = [p for block in factors for p in block]
    assert len(primes)==len(set(primes))
    assert all(p < 2**64 and p%4==1 and isprime(p) for p in primes+[k*k+1])
    norms = {e:prod(block) for e,block in zip(edges,factors)}
    assert all(norms[e]==origin+step*b for e,b in zip(edges,shifts))
    stars = []
    for i in range(5):
        ti = int(T[i].subs({P:origin, q:step}))
        star = (k*ti-R[i],ti)
        stars.append(star)
        ni = prod(norms[e] for e in incident[i])
        defect = star[0]**2+star[1]**2-ni
        assert gcd(*star)==1 and 0 < (1000*abs(defect))**5 < ni**3
    for i,j in edges:
        det = stars[i][0]*stars[j][1]-stars[i][1]*stars[j][0]
        assert det==(R[j]-R[i])*prod(norms[e] for e in incident[i]&incident[j])
    print("Numerical parent verified by exact inequalities and prime factorizations.")


def check_pair_only():
    x = (-27, -270, 135, 120, 40)
    S, e2, e3, e4, e5 = [sum(prod(t) for t in combinations(x,k)) for k in range(1,6)]
    assert e4 == 0 and isqrt(-S*e5)**2 == -S*e5
    b = {e:x[e[0]]*x[e[1]]*(S-x[e[0]]-x[e[1]]) for e in edges}
    assert len(set(b.values()) | {0})==11
    h = [t**3-S*t*t+e2*t for t in x]
    L, B = isqrt(-S*e5)*q, Rational(3,2)*e5*q*q
    Y = [P+q*t for t in h]
    X = [P*y+B*t for y,t in zip(Y,x)]
    for i,j in edges:
        assert expand(X[i]*L*Y[j]-X[j]*L*Y[i]-L*B*(x[i]-x[j])*(P+q*b[i,j])) == 0
    leading = []
    for i in range(5):
        n = prod(P+q*b[e] for e in edges if i in e)
        dx = Poly(expand(X[i]**2+(L*Y[i])**2-n),P,q)
        assert all(a<=1 and a+b==4 for (a,b),_ in dx.terms())
        c = e5*x[i]*(e2*x[i]-S*e2-e3)
        assert dx.coeff_monomial(P*q**3) == c
        leading.append(c)
    assert min(leading)<0<max(leading)
    print("Uncomplemented: ten determinants; five shell identities; both defect signs.")

if __name__ == "__main__":
    check_balanced()
    check_pair_only()
