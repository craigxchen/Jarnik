"""Exact jet lattice/image indices and full-profile product exponents."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd
from random import Random

from check_oriented_contact_lattice_rigidity import (
    add, bezout, conj, exact, lin, mul, norm, power,
)


def egcd(a, b):
    r, r1, s, s1, t, t1 = a, b, 1, 0, 0, 1
    while r1:
        q = r//r1
        r, r1, s, s1, t, t1 = r1, r-q*r1, s1, s-q*s1, t1, t-q*t1
    if r < 0:
        r, s, t = -r, -s, -t
    return r, s, t


def add_congruence(B, row, modulus):
    n = len(B)
    v = [sum(row[i]*B[i][j] for i in range(n)) for j in range(n)]
    for j in range(1, n):
        a, b = v[0], v[j]
        if b == 0:
            continue
        g, s, t = egcd(a, b)
        old0, oldj = [B[i][0] for i in range(n)], [B[i][j] for i in range(n)]
        for i in range(n):
            B[i][0] = s*old0[i]+t*oldj[i]
            B[i][j] = (-b//g)*old0[i]+(a//g)*oldj[i]
        v[0], v[j] = g, 0
    factor = modulus//gcd(modulus, v[0])
    for i in range(n):
        B[i][0] *= factor
    return factor


def transform_matrix(d, lam):
    r, s = lam
    u, v = bezout(r, s)
    # Columns (r,s),(-v,u), with determinant one, send L_lambda to Y.
    M = [[0]*(d+1) for _ in range(d+1)]
    for j in range(d+1):
        for a in range(d-j+1):
            for b in range(j+1):
                M[a+b][j] += (comb(d-j,a)*r**(d-j-a)*(-v)**a
                              *comb(j,b)*s**(j-b)*u**b)
    return M


def evaluate(P, z, w):
    d, value = len(P)-1, (0, 0)
    for k, c in enumerate(P):
        term = mul(power(w,d-k), power((-z[0],-z[1]),k))
        value = add(value, lin(c,term,0,(0,0)))
    return value


def jet_indices():
    z, w = (1,0), (2,1)
    specs = [((2,1),2), ((3,2),1)]
    contacts, delta = [], (1,0)
    for pi, e in specs:
        g = power(pi,e)
        contacts.append((norm(g),(g[0]-2*g[1],g[1])))
        delta = mul(delta,g)
    D = norm(delta)
    count = 0
    for d in range(1,7):
        for h in range(1,d+1):
            B = [[int(i==j) for j in range(d+1)] for i in range(d+1)]
            index, constraints = 1, []
            for modulus, lam in contacts:
                M = transform_matrix(d,lam)
                for k in range(h):
                    n = modulus**(h-k)
                    index *= add_congruence(B,M[k],n)
                    constraints.append((M[k],n))
            assert index == D**(h*(h+1)//2)
            values = []
            for j in range(d+1):
                P = [B[i][j] for i in range(d+1)]
                assert all(sum(a*b for a,b in zip(row,P))%n==0 for row,n in constraints)
                value = evaluate(P,z,w)
                exact(value,power(delta,h))
                values.append(value)
                size = sum(F(c*c,comb(d,k)) for k,c in enumerate(P))
                assert norm(value) <= size*(norm(z)+norm(w))**d
            image_index = 0
            for a,b in combinations(values,2):
                image_index = gcd(image_index,abs(a[0]*b[1]-a[1]*b[0]))
            assert image_index == D**h
            count += 1
    print(str(count)+' jet lattices: full prime-power coefficient indices and exact Gaussian image indices pass (d<=6, including depth two at 5).')


def clusters(m):
    out = list(range(3,m))
    result = []
    for mask in range(2**len(out)):
        T = {out[j] for j in range(len(out)) if (mask>>j)&1}
        if T:
            result.append(T)
        for a in range(3):
            result.append(T|{a})
    return result


def continuous_min(m):
    r = m-3
    candidates = [(F(0),F(0)),(F(1),F(0))]
    candidates += [(1-F(l,k),F(1,k)) for k in range(1,r+1) for l in range(k+1)]
    def objective(ab):
        a,b = ab
        return sum(comb(r,k)*(abs(k*b-1)+3*abs(a+k*b-1)) for k in range(r+1))-1
    best = min(candidates,key=objective)
    return objective(best),best


def product_exponents():
    rng = Random(581739)
    expected = {4:F(1),5:F(3),6:F(7),7:F(15),8:F(31),9:F(61),
                10:F(368,3),11:F(683,3),12:F(455)}
    for m in range(4,13):
        A, C = 2**(m-3), clusters(m)
        cm,best = continuous_min(m)
        assert cm == expected[m] and cm >= 1
        for _ in range(50):
            ns = [rng.randrange(5) for _ in range(m)]
            h,d = rng.randrange(1,6),sum(ns)
            counts = [sum(ns[i] for i in S) for S in C]
            height = A*sum(ns[3:])+2*sum(max(h-v,0) for v in counts)
            threshold = h*(4*A-1)-d*A
            gap = sum(abs(v-h) for v in counts)
            assert height-threshold == gap and gap >= h*cm
            shared = [sum(ns[i] for i in S) for S in C if len(S)>1]
            height2 = A*sum(ns[3:])+2*sum(max(h-v,0) for v in shared)
            threshold2 = h*(4*A-m-1)-d*A
            assert height2-threshold2 == d+sum(abs(v-h) for v in shared)
        a,b = best
        h = a.denominator*b.denominator
        ns = [int(h*a)]*3+[int(h*b)]*(m-3)
        assert sum(abs(sum(ns[i] for i in S)-h) for S in C) == h*cm
        r = m-3
        assert cm <= 8*comb(r-1,(r-1)//2)-1
    print('450 weighted source profiles pass exact full/private-deleted exponent identities; all nine convex minima and integer attainments pass.')


if __name__ == '__main__':
    jet_indices()
    product_exponents()
