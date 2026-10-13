"""Independent exact verification of the discrete-log collision computation (demands7.py).

Check A (all Q <= 4M, all four units): compute zhat_x mod Q by direct Gaussian multiplication of
pi_j or conj(pi_j) (no discrete logs) and compare, for every pair and unit, "zhat_x = i^s zhat_y mod Q"
with the dlog criterion kappa_x - kappa_y = s*iota mod m_Q.
Check B (big integers, M = 12 and M = 20): compute zhat_x exactly in Z[i] (hundreds of digits), and for
every pair and unit u compute the largest power q^a <= 4M dividing zhat_x - u*zhat_y (both coordinates);
compare the resulting exact collision moduli with demands7.collision_moduli (unit 1 and any unit).
Also checks Norm(zhat_x) = N exactly.
Usage: python3 verify_direct.py
"""
import random, sys
import numpy as np
from gauss7 import *
from demands7 import *

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]   # i^s, s = 0..3


def gmul_exact(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def run(q0, b, variant, trial):
    M = q0 + 1
    S = profile_construct(q0, b) if variant == 'construct' else profile_lower(q0, b)
    check_profile(S)
    Ap = (1 + S) // 2
    r = S.shape[1]
    rng = random.Random(7919 * trial + q0)
    start = rng.randrange(10 ** 7, 10 ** 9)
    ps = primes_1mod4_window(start, r)
    gA, gB = gaussian_primes(ps, rng)
    X = 4 * M
    classes = point_classes(Ap, gA, gB, ps, X)
    nQ = 0
    # ---- Check A
    for (q, a, Q, m, iota, kappa) in classes:
        z = []
        for x in range(M):
            acc = (1, 0)
            for j in range(r):
                f = (gA[j] % Q, gB[j] % Q) if Ap[x, j] else (gA[j] % Q, (-gB[j]) % Q)
                acc = gmul(acc, f, Q)
            z.append(acc)
        for x in range(M):
            for y in range(M):
                if x == y:
                    continue
                for s, u in enumerate(UNITS):
                    direct = z[x] == gmul((u[0] % Q, u[1] % Q), z[y], Q)
                    viadlog = (kappa[x] - kappa[y] - s * iota) % m == 0
                    assert direct == viadlog, (q, a, x, y, s)
        nQ += 1
    m1, many = collision_moduli(classes, X, M)
    out = "%s q0=%d b=%d trial=%d: check A passed on %d prime powers <= %d, %d ordered pairs x 4 units" % (
        variant, q0, b, trial, nQ, X, M * (M - 1))
    # ---- Check B (exact big integers)
    if M <= 20:
        N = 1
        for p in ps:
            N *= p
        zs = []
        for x in range(M):
            acc = (1, 0)
            for j in range(r):
                f = (gA[j], gB[j]) if Ap[x, j] else (gA[j], -gB[j])
                acc = gmul_exact(acc, f)
            assert acc[0] ** 2 + acc[1] ** 2 == N
            zs.append(acc)
        pp = odd_prime_powers_upto(X)
        for x in range(M):
            for y in range(M):
                if x == y:
                    continue
                mod1, modany = 1, 1
                for s, u in enumerate(UNITS):
                    uz = gmul_exact(u, zs[y])
                    d = (zs[x][0] - uz[0], zs[x][1] - uz[1])
                    assert d != (0, 0)
                    mm = 1
                    for (qq, aa, QQ) in pp:
                        if d[0] % QQ == 0 and d[1] % QQ == 0:
                            mm *= qq
                    if s == 0:
                        mod1 = mm
                    modany *= mm        # at most one unit per prime power, so the product is the any-unit modulus
                assert mod1 == m1[x, y], (x, y, mod1, m1[x, y])
                assert modany == many[x, y], (x, y, modany, many[x, y])
        out += "; check B (exact Z[i], N with %d digits) passed" % (int(N.bit_length() * 0.30103) + 1)
    print(out)
    sys.stdout.flush()


if __name__ == '__main__':
    for t in range(3):
        run(11, 8, 'construct', t)
        run(11, 33, 'lower', t)
    run(19, 8, 'construct', 0)
    run(19, 33, 'lower', 0)
    run(23, 8, 'construct', 0)
    run(31, 8, 'construct', 0)
    print("ALL DIRECT CHECKS PASSED")
