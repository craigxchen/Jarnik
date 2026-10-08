"""Exact finite audit of residual norm fibers and split-quadratic identities."""
from collections import defaultdict
from fractions import Fraction


def mul(z, w):
    return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]


def conj(z):
    return z[0], -z[1]


def add(z, w):
    return z[0]+w[0], z[1]+w[1]


def sub(z, w):
    return z[0]-w[0], z[1]-w[1]


def scale(k, z):
    return k*z[0], k*z[1]


def norm(z):
    return z[0]**2+z[1]**2


reps = defaultdict(list)
for x in range(-5, 6):
    for y in range(-5, 6):
        if 0 < x*x+y*y <= 30:
            reps[x*x+y*y].append((x,y))
coefficients = [(1,0), (1,1), (2,1)]
checks = fibers_checked = 0
for Q, betas in reps.items():
    for M, ws in reps.items():
        N = Q*M
        anchor = mul(betas[0], ws[0])
        for a in coefficients:
            for b in coefficients:
                normal = mul(a, conj(b))
                D = norm(b)*M+norm(a)*Q
                fibers = defaultdict(set)
                for beta in betas:
                    for w in ws:
                        z = mul(beta,w)
                        residual = sub(mul(b,w), mul(a,conj(beta)))
                        assert mul(beta,residual) == sub(mul(b,z),scale(Q,a))
                        assert norm(residual) == D-2*mul(conj(normal),z)[0]
                        assert (norm(residual)-D) % 2 == 0
                        fibers[norm(residual)].add(z)
                        u = sub(scale(Q,w),mul(anchor,conj(beta)))
                        v = add(scale(Q,w),mul(anchor,conj(beta)))
                        assert mul(u,conj(v))[0] == 0
                        assert norm(u)+norm(v) == 4*Q*N
                        assert mul(u,add(z,anchor)) == mul(v,sub(z,anchor))
                        denominator = add(z, anchor)
                        if denominator == (0, 0):
                            # The omitted projective tangent-half-angle point.
                            assert v == (0, 0)
                            assert u != (0, 0)
                        else:
                            assert v != (0, 0)
                            residual_ratio = mul(u, conj(v))
                            anchor_ratio = mul(sub(z, anchor), conj(denominator))
                            assert residual_ratio[0] == 0
                            assert anchor_ratio[0] == 0
                            assert Fraction(residual_ratio[1], norm(v)) == \
                                Fraction(anchor_ratio[1], norm(denominator))
                        checks += 1
                for residual_norm, products in fibers.items():
                    assert len(products) <= (1 if residual_norm == 0 else 2)
                    fibers_checked += 1
print(f"Verified {checks} exact tuples and {fibers_checked} residual-norm fibers.")
