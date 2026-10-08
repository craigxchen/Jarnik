"""Exact full-incidence span/saturation and quotient-Jacobian checks."""
from math import comb, gcd, lcm
from random import Random

from check_higher_degree_gram_jet_barrier import add_congruence, egcd, transform_matrix, clusters
from check_oriented_contact_lattice_rigidity import crt, add, mul, norm, power


def compositions(total, count):
    if count == 1:
        yield (total,)
        return
    for n in range(total+1):
        for rest in compositions(total-n, count-1):
            yield (n,)+rest


def multiply(P, Q):
    out = [0]*(len(P)+len(Q)-1)
    for i,a in enumerate(P):
        for j,b in enumerate(Q):
            out[i+j] += a*b
    return out


def column_basis(A):
    A = [row[:] for row in A]
    r, n = len(A), len(A[0])
    for k in range(r):
        j = next(j for j in range(k,n) if A[k][j])
        for row in A:
            row[k], row[j] = row[j], row[k]
        for j in range(k+1,n):
            a,b = A[k][k],A[k][j]
            if not b:
                continue
            g,s,t = egcd(a,b)
            left,right = [row[k] for row in A],[row[j] for row in A]
            for i in range(r):
                A[i][k] = s*left[i]+t*right[i]
                A[i][j] = (-b//g)*left[i]+(a//g)*right[i]
        if A[k][k] < 0:
            for row in A:
                row[k] *= -1
    assert all(A[i][j] == 0 for i in range(r) for j in range(r,n))
    B = [row[:r] for row in A]
    for k in range(r-1,-1,-1):
        for j in range(k):
            q = B[k][j]//B[k][k]
            for i in range(r):
                B[i][j] -= q*B[i][k]
    return B


def membership(B, vector):
    v = vector[:]
    for k in range(len(B)):
        if v[k] % B[k][k]:
            return False
        c = v[k]//B[k][k]
        for i in range(k,len(B)):
            v[i] -= c*B[i][k]
    return all(t == 0 for t in v)


def source_data(depth):
    C = clusters(5)
    primes = [5,13,17,29,37,41,53,61,73,89,97,101,109,113,137]
    moduli = [p**(depth if j==0 else 1) for j,p in enumerate(primes)]
    r = [0,5,10]
    for i in [3,4]:
        residues = []
        for j,S in enumerate(C):
            if i in S:
                anchor = next((a for a in S if a < 3),None)
                value = r[anchor] if anchor is not None else 0
            else:
                value = 5 if j == 0 else i+3
            residues.append(value)
        r.append(crt(residues,moduli))
    for S,n in zip(C,moduli):
        assert all((r[i]-r[min(S)])%n == 0 for i in S)
    return r,C,moduli


def span_checks():
    cases = 0
    for depth in [2,5]:
        r,C,moduli = source_data(depth)
        delta = lcm(abs(r[0]-r[1]),abs(r[0]-r[2]),abs(r[1]-r[2]))
        assert delta == 10
        D = 1
        for n in moduli:
            D *= n
        for d,h in [(1,1),(2,1),(2,2),(3,2)]:
            generators = []
            for ns in compositions(d,5):
                P,M = [1],1
                for i,n in enumerate(ns):
                    for _ in range(n):
                        P = multiply(P,[-1,r[i]])
                for S,norm_e in zip(C,moduli):
                    M *= norm_e**max(h-sum(ns[i] for i in S),0)
                generators.append([M*c for c in P])
            source_B = column_basis([list(row) for row in zip(*generators)])
            source_index = 1
            for i in range(d+1):
                source_index *= source_B[i][i]
            jet_B = [[int(i==j) for j in range(d+1)] for i in range(d+1)]
            jet_index,conditions = 1,[]
            for S,norm_e in zip(C,moduli):
                lam = (r[min(S)],1)
                T = transform_matrix(d,lam)
                for k in range(h):
                    modulus = norm_e**(h-k)
                    jet_index *= add_congruence(jet_B,T[k],modulus)
                    conditions.append((T[k],modulus))
            assert jet_index == D**(h*(h+1)//2)
            assert source_index % jet_index == 0
            ratio = source_index//jet_index
            assert delta**(d*(d+1)) % ratio == 0
            for P in generators:
                assert all(sum(a*b for a,b in zip(row,P))%n == 0 for row,n in conditions)
            for j in range(d+1):
                target = [delta**d*jet_B[i][j] for i in range(d+1)]
                assert membership(source_B,target)
            print('depth',depth,'degree/order',(d,h),'exact saturation index',ratio)
            cases += 1
    print(str(cases)+' full fifteen-cluster spans satisfy delta^d Gamma subset Lambda, including shallow outside collisions at a deep contact prime.')


def jacobian_checks():
    rng = Random(85822)
    for _ in range(50):
        n = rng.randrange(-20,21)
        z,w = (1,0),(n,1)
        H2 = norm(z)+norm(w)
        for d in range(1,9):
            rr,ri,ii = 0,0,0
            for k in range(d+1):
                value = mul(power(w,d-k),power((-z[0],-z[1]),k))
                weight = comb(d,k)
                rr += weight*value[0]**2
                ri += weight*value[0]*value[1]
                ii += weight*value[1]**2
            assert 4*(rr*ii-ri*ri) == H2**(2*d)-(H2*H2-4)**d
            assert rr+ii == H2**d
    print('400 exact Bombieri evaluation Gram determinants pass the quotient-Jacobian formula.')


if __name__ == '__main__':
    span_checks()
    jacobian_checks()
