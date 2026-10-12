"""Exact check of the two-logarithm identity (Section 2 of flipped.md).

For every nonconstant label a:   prod_x zeta_x^{H(x,a)} = +- u_a^M v_a^{-2},
   zeta_x = z_x/conj(z_x), u_a = G_a/conj(G_a), v_a = P_a/conj(P_a),
   G_a = prod_{a_j=a} pi_j^(sigma_j),  P_a = prod_j pi_j^(sigma_j c_j(a)),
   c_j(a) = sum_{x in F_j} H(x,a) H(x,a_j).
Equivalently A_a * conj(Gamma_a) is real or purely imaginary, where
   A_a = prod_{H(x,a)=1} z_x * prod_{H(x,a)=-1} conj(z_x),   Gamma_a = G_a^M conj(P_a)^2.
Pair form: (A_a conj(A_b)) * conj(Gamma_ab) real or imaginary, Gamma_ab = (G_a conj(G_b))^M conj(Q_ab)^4,
   Q_ab = prod_j pi_j^(sigma_j (c_j(a)-c_j(b))/2).
Also: the conj-primitive parts of Gamma_a, Gamma_ab are non-units off the axes and diagonals (so the
linear forms are nonzero), and the naive height of gamma/conj(gamma) is max(N, 2|Re gamma^2|) <= 2N.
All arithmetic exact (Python integers).
"""
import math, random, sys
from fcommon import *


def label_data(P, a):
    H, M = P.H, P.M
    Ga = gprod([gsignpow(c["pi"], c["sigma"]) for c in P.cols if c["label"] == a])
    cj = [sum(H[x][a] * H[x][c["label"]] for x in c["F"]) for c in P.cols]
    Pa = gprod([gsignpow(c["pi"], c["sigma"] * cj[j]) for j, c in enumerate(P.cols)])
    return Ga, Pa, cj


def real_or_imag(z):
    return z[0] == 0 or z[1] == 0


def check_profile(P, rng, npairs=12):
    H, M = P.H, P.M
    Z = [P.point(x) for x in range(M)]
    n = None
    for z in Z:
        if n is None:
            n = gnorm(z)
        assert gnorm(z) == n
    ok = 0
    data = {}
    for a in range(1, M):
        Ga, Pa, cj = label_data(P, a)
        data[a] = (Ga, Pa, cj)
        A = gprod([Z[x] if H[x][a] == 1 else gconj(Z[x]) for x in range(M)])
        Gam = gmul(gpow(Ga, M), gpow(gconj(Pa), 2))
        assert real_or_imag(gmul(A, gconj(Gam))), ("identity fails", a)
        # nonvanishing of the form: conj-primitive part not on axes / diagonals
        exps = [c["sigma"] * ((M if c["label"] == a else 0) - 2 * cj[j]) for j, c in enumerate(P.cols)]
        g = conj_primitive_part(exps, [c["pi"] for c in P.cols])
        assert g[0] != 0 and g[1] != 0 and abs(g[0]) != abs(g[1]), "form vanishes?"
        ok += 1
    labels = list(range(1, M))
    for _ in range(npairs):
        a, b = rng.sample(labels, 2)
        Ga, Pa, ca = data[a]
        Gb, Pb, cb = data[b]
        diff = [ca[j] - cb[j] for j in range(len(P.cols))]
        assert all(d % 2 == 0 for d in diff)
        Q = gprod([gsignpow(c["pi"], c["sigma"] * diff[j] // 2) for j, c in enumerate(P.cols)])
        U = gmul(Ga, gconj(Gb))
        Gam = gmul(gpow(U, M), gpow(gconj(Q), 4))
        Aa = gprod([Z[x] if H[x][a] == 1 else gconj(Z[x]) for x in range(M)])
        Ab = gprod([Z[x] if H[x][b] == 1 else gconj(Z[x]) for x in range(M)])
        lhs = gmul(gmul(Aa, gconj(Ab)), gconj(Gam))
        assert real_or_imag(lhs), ("pair identity fails", a, b)
        # Q's norm = prod over flipped rows in D_ab of their flipped primes (one flip per column)
        Dab = {x for x in range(M) if H[x][a] != H[x][b]}
        expect = 1
        for c in P.cols:
            for x in c["F"]:
                if x in Dab:
                    expect *= c["p"]
        if all(len(c["F"]) <= 1 for c in P.cols):
            assert gnorm(Q) == expect
        ok += 1
    return ok


def naive_height_check(rng, trials=2000):
    cnt = 0
    for _ in range(trials):
        k = rng.randint(1, 4)
        ps = rng.sample(split_primes_from(5, 40), k)
        g = gprod([gsignpow(gauss_prime(p), rng.choice([1, -1]) * rng.randint(1, 3)) for p in ps])
        N = gnorm(g)
        re2 = gmul(g, g)[0]
        assert math.gcd(N, 2 * re2) == 1
        assert max(N, abs(2 * re2)) <= 2 * N
        cnt += 1
    return cnt


def main():
    rng = random.Random(20261009)
    total = 0
    cases = [("Walsh-8", sylvester(3)), ("Paley-12", paley_I(11)), ("Walsh-16", sylvester(4)),
             ("Paley-20", paley_I(19))]
    for name, H in cases:
        M = len(H)
        for trial in range(4):
            primes = split_primes_from(rng.randint(5, 400), 5 * (M - 1))
            rng.shuffle(primes)
            # random capacity-two assignment: a surjection-ish map rows -> labels with each label <= 2
            labels = list(range(1, M)) + [rng.randrange(1, M)]
            rng.shuffle(labels)
            P, _ = fully_flipped_profile(H, 5, primes, assignment=labels[:M], rng=rng)
            total += check_profile(P, rng)
        # general modifications: several flips per column, flips on content columns
        for trial in range(2):
            primes = split_primes_from(rng.randint(5, 400), 5 * (M - 1) + 3)
            P, _ = fully_flipped_profile(H, 5, primes, rng=rng)
            for c in P.cols:
                if rng.random() < 0.15:
                    c["F"].add(rng.randrange(M))
            for extra in primes[5 * (M - 1):]:
                P.cols.append({"label": 0, "sigma": 1, "p": extra, "pi": gauss_prime(extra),
                               "F": {rng.randrange(M)}})
            total += check_profile(P, rng)
        print(f"{name}: identities verified")
    total += naive_height_check(rng)
    # negative control: wrong exponents must fail (the identity is not vacuous)
    H = paley_I(11); M = 12
    P, _ = fully_flipped_profile(H, 5, split_primes_from(50, 55), rng=rng)
    Z = [P.point(x) for x in range(M)]
    for a in range(1, M):
        Ga, Pa, cj = label_data(P, a)
        A = gprod([Z[x] if H[x][a] == 1 else gconj(Z[x]) for x in range(M)])
        for e1, e2 in [(M - 2, 2), (M, 1), (M, 0), (M, 4)]:
            Gam = gmul(gpow(Ga, e1), gpow(gconj(Pa), e2))
            assert not real_or_imag(gmul(A, gconj(Gam))), "negative control unexpectedly passed"
    print("negative control: wrong exponents fail as they should")
    print("TOTAL exact checks passed:", total)
    print("PASS")


if __name__ == "__main__":
    main()
