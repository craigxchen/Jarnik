"""Max number of lattice points of norm q (squarefree product of split primes) on an arc of
length C q^{1/4} of the circle x^2+y^2=q.  Own code; exact Gaussian products, float angles
(the angular windows are ~q^{-1/4}, far above double precision), plus an exact verification
of the Gaussian triangle identity z_j - z_i = 2 i eta H_ij t_ij on every close pair/triple."""
import sys, random, itertools, cmath
from math import isqrt, prod, log, atan2, pi
from sympy import primerange


def two_squares(p):
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        b = isqrt(b2)
        if b * b == b2:
            return (a, b)


def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def points(primes):
    gens = [two_squares(p) for p in primes]
    pts = []
    for choice in itertools.product((0, 1), repeat=len(primes)):
        z = (1, 0)
        for g, c in zip(gens, choice):
            z = gmul(z, g if c == 0 else (g[0], -g[1]))
        for u in ((1, 0), (0, 1), (-1, 0), (0, -1)):
            pts.append((gmul(z, u), choice, u))
    return pts


def max_arc(primes, C):
    q = prod(primes)
    pts = points(primes)
    for (z, _, _) in pts:
        assert z[0] ** 2 + z[1] ** 2 == q
    R = q ** 0.5
    width = C * q ** 0.25 / R  # angular width of the arc
    ang = sorted((atan2(z[1], z[0]) % (2 * pi), k) for k, (z, _, _) in enumerate(pts))
    n = len(ang)
    best, bestset = 1, None
    j = 0
    ext = ang + [(a + 2 * pi, k) for a, k in ang]
    for i in range(n):
        while j < i + n and ext[j][0] - ext[i][0] <= width:
            j += 1
        if j - i > best:
            best = j - i
            bestset = [pts[ext[t][1]] for t in range(i, j)]
    return best, bestset, q


def check_pair_identity(primes, P1, P2):
    """same unit eta: z2 - z1 = 2 i eta H t with t a nonzero rational integer."""
    (z1, c1, u1), (z2, c2, u2) = P1, P2
    if u1 != u2:
        return None
    gens = [two_squares(p) for p in primes]
    H = (1, 0)
    for g, a, b in zip(gens, c1, c2):
        if a == b:
            H = gmul(H, g if a == 0 else (g[0], -g[1]))
    d = (z2[0] - z1[0], z2[1] - z1[1])
    # divide d by 2 i u1 H : multiply by conj and divide by norms
    den = gmul(gmul((0, 2), u1), H)
    nd = den[0] ** 2 + den[1] ** 2
    num = gmul(d, (den[0], -den[1]))
    assert num[0] % nd == 0 and num[1] % nd == 0
    t = (num[0] // nd, num[1] // nd)
    assert t[1] == 0 and t[0] != 0, t
    return t[0]


if __name__ == '__main__':
    seed, s, P, C, trials = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
    rng = random.Random(seed)
    pool = [p for p in primerange(P, 2 * P) if p % 4 == 1]
    dist = {}
    worst = None
    npairs = 0
    for tr in range(trials):
        primes = sorted(rng.sample(pool, s))
        M, S, q = max_arc(primes, C)
        dist[M] = dist.get(M, 0) + 1
        if S:
            for A, B in itertools.combinations(S, 2):
                if check_pair_identity(primes, A, B) is not None:
                    npairs += 1
        if worst is None or M > worst[0]:
            worst = (M, primes)
    W = s * log(1.5 * P)
    print(f'circle s={s} P={P} K={len(pool)} C={C} trials={trials} dist={dict(sorted(dist.items()))} '
          f'worst={worst} same-unit pairs checked={npairs} logq/loglogq={W/log(W):.2f}')
