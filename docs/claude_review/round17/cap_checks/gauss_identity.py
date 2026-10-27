# Verify (exactly, in Z[zeta] via integer polynomial arithmetic mod Phi_q) the slab-coordinate identities
#   (1 - i sqrt q)/2 * hat k(1) = sum_a zeta^a S_a ,  S_a = k(Pp_a)   (Pp_a = {a} u (a+QR))
#   (1 + i sqrt q)/2 * hat k(1) = sum_a zeta^a T_a ,  T_a = k(P_a)    (P_a  = {a} u (a+NQR))
# using i sqrt q = g = sum_x chi(x) zeta^x (Gauss sum, q = 3 mod 4).  Equivalent integer identities in Z[zeta]:
#   (1 - g) * hat k(1) = 2 * sum_a zeta^a S_a ,   (1 + g) * hat k(1) = 2 * sum_a zeta^a T_a .
# Also: S_a + T_a = K + k_a.  Random integer k (not only ternary); identities are linear so this is a proof
# check of the algebra for each q (the general proof is in cap.md §2).
import random, sys
from sdclib import chi, patterns
def mulmod(f, g, q):
    h = [0] * q
    for i, a in enumerate(f):
        if a:
            for j2, b in enumerate(g):
                if b: h[(i + j2) % q] += a * b
    return h
def canon(f, q):
    # reduce mod 1 + z + ... + z^{q-1}: subtract f[q-1] * (all ones) to make top coeff 0
    c = f[q - 1]
    return [a - c for a in f]
out = []
for q in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 103, 107]:
    P, Pp = patterns(q)
    g = [chi(x, q) for x in range(q)]
    one = [1] + [0] * (q - 1)
    ok = True
    for trial in range(20):
        k = [random.randint(-5, 5) for _ in range(q)]
        hk = k[:]  # hat k(1) = sum k_x zeta^x
        S = [sum(k[x] for x in Pp[a]) for a in range(q)]
        T = [sum(k[x] for x in P[a]) for a in range(q)]
        K = sum(k)
        ok &= all(S[a] + T[a] == K + k[a] for a in range(q))
        lhs1 = canon(mulmod([o - gg for o, gg in zip(one, g)], hk, q), q)
        rhs1 = canon([2 * s for s in S], q)
        lhs2 = canon(mulmod([o + gg for o, gg in zip(one, g)], hk, q), q)
        rhs2 = canon([2 * t for t in T], q)
        ok &= (lhs1 == rhs1) and (lhs2 == rhs2)
    # g^2 = -q check
    g2 = canon(mulmod(g, g, q), q)
    ok &= (g2 == canon([-q] + [0] * (q - 1), q))
    line = f"q={q}: identities (1-g)hat k(1)=2 sum zeta^a S_a, (1+g)hat k(1)=2 sum zeta^a T_a, S+T=K+k, g^2=-q : {ok}"
    print(line, flush=True); out.append(line)
open('gauss_identity_out.txt', 'w').write('\n'.join(out) + '\n')
