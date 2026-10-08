"""Exact check of Theorem B(i): for common-unit lattice points z_a on x^2+y^2=N and anchored
numerators P_i (P_0=1), z_i - z_j = -2i z_0 det(P_i,P_j)/(conj(P_i) conj(P_j))."""
import random, sys
sys.path.insert(0, '.')
from gi import *
random.seed(3)
PR = split_primes(80); PI = {p: gaussian_prime_above(p) for p in PR}
ok = 0
for trial in range(200):
    ps = random.sample(PR, 3); es = [random.randint(1, 3) for _ in ps]
    allocs = [tuple(random.randint(0, e) for e in es) for _ in range(5)]
    if len(set(allocs)) < 5:
        continue
    def pt(a):
        z = (1, 0)
        for p, ai, e in zip(ps, a, es):
            z = mul(z, mul(gpow(PI[p], ai), gpow(conj(PI[p]), e - ai)))
        return z
    z = [pt(a) for a in allocs]
    P = [(1, 0)]
    for a in allocs[1:]:
        q = (1, 0)
        for idx, p in enumerate(ps):
            d = a[idx] - allocs[0][idx]
            q = mul(q, gpow(PI[p], d) if d > 0 else gpow(conj(PI[p]), -d))
        P.append(q)
    for i in range(5):
        for j in range(5):
            if i == j:
                continue
            det = P[i][0] * P[j][1] - P[j][0] * P[i][1]
            lhs = mul(sub(z[i], z[j]), mul(conj(P[i]), conj(P[j])))
            rhs = mul((0, -2), smul(det, z[0]))
            assert lhs == rhs
            ok += 1
print("Theorem B(i) verified exactly on", ok, "ordered pairs")
