"""Independent re-verification of a joint witness pickle (round 14).  Checks, in exact integers:
(1) parity of every level-1 label and of its level-2 lift; lift reduces to the level-1 label;
(2) no level-2 collision of any v in V_K at the non-G moduli (so no deeper collisions);
(3) every budget m_(v,s)^2 <= Norm(Z_v) (2 if s odd), with genuine all-level demand at q in G;
(4) the label at q0 lies outside every per-column-conjugation class {u (eps o L) mod n}.
Also recomputes the forced set and checks the witness covers exactly F_K^ex."""
import sys, pickle, math, itertools
from rlib14 import *
from forced import forced_set
w = pickle.load(open(sys.argv[1], "rb"))
K, ps, G, q0, labels, lev2 = w['K'], w['ps'], w['G'], w['q0'], w['labels'], w['lev2']
r = len(ps); pis = gaussian_primes(ps); F = forced_set(ps, K, 60000)
assert sorted(list(labels) + list(G)) == sorted(F), "moduli mismatch"
bad = 0
for q in labels:
    m = q - chi4(q); nv = nu(ps, q)
    assert all(labels[q][j] % 2 == nv[j] for j in range(r)), ("parity", q)
    assert all(lev2[q][j] % m == labels[q][j] % m for j in range(r)), ("lift", q)
    M2 = m * q
    for v in VK(r, K):
        t2 = sum(a * b for a, b in zip(v, lev2[q])) % M2
        assert t2 % (M2 // 4) != 0, ("level-2 collision", q, v)
viol = 0; checked = 0
for v in VK(r, K):
    z = Zv(pis, v); nz = z[0] ** 2 + z[1] ** 2
    for s in range(4):
        mm = 1
        for q in G: mm *= q ** vq(Ds(z, s), q)
        for q, ell in labels.items():
            m = q - chi4(q); n = m // 4
            if sum(a * b for a, b in zip(v, ell)) % m == s * n: mm *= q
        checked += 1
        if mm * mm > nz * (2 if s % 2 else 1): viol += 1
m0 = q0 - chi4(q0); n0 = m0 // 4; _, L0 = genuine_labels(ps, q0)
c0 = tuple(x % n0 for x in labels[q0])
inany = any(tuple((u * e * x) % n0 for e, x in zip(eps, L0)) == c0
            for eps in itertools.product([1, -1], repeat=r) for u in range(1, n0) if math.gcd(u, n0) == 1)
print(f"{sys.argv[1]}: K={K} G={G} q0={q0}: budgets checked={checked} violations={viol}; q0 label in some eps-class: {inany}")
