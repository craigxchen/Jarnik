"""Exact bounded checks for the compatible Paley coordinate-gap relaxation."""

from fractions import Fraction as F
from itertools import product
from random import Random

from check_general_four_row_quartic_slack_obstruction import paley


def profile(q):
    h, _ = paley(q)
    m, r = q+1, 5*q
    labels = [1+j%q for j in range(r)]
    s = [[h[i][label] for label in labels] for i in range(m)]
    flipped = [q]+list(range(q))
    for i,j in enumerate(flipped):
        s[i][j] *= -1
    return h,s,labels,flipped


def budget(c,h,s,labels,flipped):
    # Compute the same physical sum using the unflipped transform plus
    # the exact correction at the one physical column assigned to each row.
    m = len(c)
    t = [sum(c[i]*h[i][a] for i in range(m)) for a in range(m)]
    assert all(x%2 == 0 for x in t)
    b = 5*sum(map(abs,t))//2
    for i,j in enumerate(flipped):
        a = labels[j]
        b += abs(t[a]//2-c[i]*h[i][a])-abs(t[a]//2)
    return b


def support_tiers():
    hmat,s,labels,flipped = profile(11)
    m,r = len(s),len(s[0])
    count = exceptional = 0
    for c in product((-1,0,1),repeat=m):
        if not any(c) or sum(c):
            continue
        h = sum(x!=0 for x in c)
        b = budget(c,hmat,s,labels,flipped)
        if h == 2:
            if c[0]:
                assert 2*b >= 5*m-4
                exceptional += 1
            else:
                assert 2*b >= 5*m
        else:
            assert h>=4 and 2*b >= 5*m+2*h-8
        assert 2*b-r >= 1
        count += 1
    assert exceptional == 2*(m-1) <= 2*m
    rng = Random(39254)
    for _ in range(1000):
        c = [rng.randrange(-6,7) for _ in range(m-1)]
        c.append(-sum(c))
        a = max(map(abs,c))
        if a<2:
            continue
        b = budget(c,hmat,s,labels,flipped)
        direct = sum(abs(sum(c[i]*s[i][j] for i in range(m)))//2 for j in range(r))
        assert b == direct and 2*b >= 3*m*a
    print(f"M={m}: {count} ternary tiers, {exceptional} exceptional pair candidates")


def rational_union_bounds():
    for m in (12,20,24,44,60,84,100):
        tail = F(1,m**5)*F(5,m)**m/(1-F(5,m**3)**m)
        q = 2+F(1,m**3)+F(16,m-2)+tail
        assert q<4
        assert 56*q/F(1000)<F(224,1000)<1
        for h in range(4,m+1):
            assert (2*m)**h*F(1,m**(2*h-3)) == F(2**h,m**(h-3))
        for a in range(2,20):
            assert 2*a+1 <= 5**(a-1)
            bound = F(5**(m*(a-1)),m**((3*a-5)*m+5))
            assert bound == F(1,m**5)*F(5,m)**m*F(5,m**3)**(m*(a-2))
    assert 2+F(1,1728)+F(8,5)+F(2,12**5)<4
    # Both sine and cosine are excluded: the final constant is 56, not 28.
    assert 2*14*2 == 56
    print("PASS: exact support/amplitude series and doubled coordinate-gap budget")


def phase_lift():
    for q in (11,19,43):
        h,s,labels,flipped = profile(q)
        m,r = len(s),len(s[0])
        phi = [F(i+1,10000*m) for i in range(m)]
        theta = [F(0) for _ in range(r)]
        for i,j in enumerate(flipped):
            a = labels[j]
            unflipped = 2*q+(a-1)
            assert all(s[k][j]-s[k][unflipped]
                       == (-2*h[i][a] if k==i else 0) for k in range(m))
            coefficient = phi[i]/(-2*h[i][a])
            theta[j] += coefficient
            theta[unflipped] -= coefficient
        assert [sum(s[i][j]*theta[j] for j in range(r)) for i in range(m)] == phi
        assert any(x==0 for x in theta)
        for i in range(1,m):
            c = [-1]+[int(k==i) for k in range(1,m)]
            v = [sum(c[k]*s[k][j] for k in range(m))//2 for j in range(r)]
            assert sum(v[j]*theta[j] for j in range(r)) == (phi[i]-phi[0])/2
    print("PASS: exact rational full-rank lift and formal character phase identities")


if __name__ == "__main__":
    support_tiers()
    rational_union_bounds()
    phase_lift()
    print("PASS: compatible-phase relaxation checks; no Gaussian integrality claimed")
