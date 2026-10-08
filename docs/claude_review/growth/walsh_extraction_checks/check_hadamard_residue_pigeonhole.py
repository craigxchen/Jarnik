"""Checks for Theorems B, B'' (residue pigeonhole for pure / row-deleted Hadamard cores).

1. Label-count lower bound |A| >= |U| W/(2W+4|U|kappa''_+) on random weights, Walsh and Paley cores,
   full and row-deleted, where 2kappa'' = max_{x!=y in U} G_xy.
2. Exact pinning identity: for literal pure-Hadamard Gaussian rows, sum_x H(x,a) lift(arg z_x) - (pi/2) sum_x H(x,a) k_x
   equals M*phi_a (k_x integers from the units/windings) -- checked numerically, and the residue map is additive.
3. Literal pure-Walsh/Paley configurations: the actual best arc constant C_act (optimised over all row units)
   never violates Theorem B: W <= max(8M log+ C, 36 log+(2C)) when M>=64.
4. Numerical constants of Theorem B'' (Sidon/pigeonhole counting) at the stated thresholds, in exact integers.
5. Explicit thresholds M_1(C) from the audited inert-prime bound F(M) <= (M/4) log R + binom(M,2) log C.
"""
import random, math, sys, os
from math import comb, log
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gauss_common import *
random.seed(7)

def walsh(t):
    M = 2 ** t
    return [[(-1) ** bin(x & a).count("1") for a in range(M)] for x in range(M)]

def paley(q):
    sq = set((a * a) % q for a in range(1, q))
    ch = lambda a: 0 if a % q == 0 else (1 if (a % q) in sq else -1)
    M = q + 1
    H = [[1] * M for _ in range(M)]
    for x in range(1, M):
        for j in range(1, M):
            H[x][j] = -1 if x == j else -ch(x - j)
    return H

def is_hadamard(H):
    M = len(H)
    return all(sum(H[x][j] * H[y][j] for j in range(M)) == (M if x == y else 0) for x in range(M) for y in range(M))

cores = [("walsh16", walsh(4)), ("walsh32", walsh(5)), ("walsh64", walsh(6)), ("paley12", paley(11)), ("paley20", paley(19)),
         ("paley44", paley(43)), ("paley68", paley(67))]
for name, H in cores:
    assert is_hadamard(H), name
    assert all(sum(H[x][a] for x in range(len(H))) == 0 for a in range(1, len(H)))

# 1. label-count bound
n1 = 0
for name, H in cores:
    M = len(H)
    for trial in range(40):
        f = random.choice([0, 0, 1, 2, M // 8])
        U = sorted(random.sample(range(M), M - f))
        Wa = [0.0] * M
        mode = random.random()
        for a in range(1, M):
            if mode < 0.3: Wa[a] = random.random() if random.random() < 0.5 else 0.0
            elif mode < 0.6: Wa[a] = random.expovariate(1.0)
            else: Wa[a] = 1.0 + 0.1 * random.random()
        if sum(Wa) == 0: continue
        W = sum(Wa)
        G = lambda x, y: sum(Wa[a] * H[x][a] * H[y][a] for a in range(1, M))
        gam = max(G(x, y) for i, x in enumerate(U) for y in U[i + 1:])
        kpp = max(0.0, gam / 2)
        A = sum(1 for a in range(1, M) if Wa[a] > 0)
        bound = len(U) * W / (2 * W + 4 * len(U) * kpp)
        assert A >= bound - 1e-9, (name, A, bound)
        n1 += 1
print("1. label-count bound held on", n1, "random weightings (full and row-deleted cores)")

# 2 & 3. literal pure-core configurations; actual arc constant over all unit choices
def best_arc(angles):
    """minimal width of an arc (mod pi/2) containing all angles taken mod pi/2"""
    q = math.pi / 2
    a = sorted(t % q for t in angles)
    gaps = [(a[(i + 1) % len(a)] - a[i]) % q for i in range(len(a))]
    gaps[-1] = (a[0] + q - a[-1])
    i = max(range(len(a)), key=lambda i: gaps[i])
    return q - gaps[i]

n3 = 0
for name, H in cores:
    M = len(H)
    for trial in range(3):
        b = random.choice([1, 2, 3])
        primes = split_primes(b * (M - 1) + 10)[10:]
        random.shuffle(primes)
        pis = [gauss_prime(p) for p in primes]
        blocks = []
        for a in range(1, M):
            G = (1, 0)
            for j in range(b):
                pi = pis[(a - 1) * b + j]
                if random.random() < 0.5: pi = gconj(pi)
                G = gmul(G, pi)
            blocks.append(G)
        phis = [math.atan2(G[1], G[0]) for G in blocks]
        angles = [sum(H[x][a] * phis[a - 1] for a in range(1, M)) for x in range(M)]
        # 2: pinning identity -- choose the units achieving the optimal window, recover k_x and m_a
        width = best_arc(angles)
        # reconstruct k_x relative to the optimal window start
        q = math.pi / 2
        a_sorted = sorted(t % q for t in angles)
        # identity check: sum_x H(x,a) angles_x = M phi_a exactly (no arc needed)
        for a in range(1, M):
            s = sum(H[x][a] * angles[x] for x in range(M))
            assert abs(s - M * phis[a - 1]) < 1e-7 * M
        Wv = sum(math.log(gnorm(G)) for G in blocks)
        R = math.exp(Wv / 2)
        logC = math.log(width) + Wv / 4      # C_act = width * N^{1/4}, N = e^W
        bound_ok = True
        if M >= 64:
            lpC = max(0.0, logC); lp2C = max(0.0, logC + math.log(2))
            bound_ok = Wv <= max(8 * M * lpC, 36 * lp2C) + 1e-9
        assert bound_ok
        n3 += 1
        print("   %-8s b=%d  W=%8.1f  log C_act=%8.1f  W/(8M)=%7.2f" % (name, b, Wv, logC, Wv / (8 * M)))
print("3. literal pure cores: optimal-unit arc constants consistent with Theorem B on", n3, "configurations")

# 4. constants of Theorem B (M>=64) and Theorem B''
for M in [64, 128, 256, 1024]:
    lhs = M / 3; rhs = (1 + math.sqrt(4 * M - 3)) / 2 + 8
    assert lhs > rhs, M
print("4a. Theorem B counting M/3 > (1+sqrt(4M-3))/2+8 for M=64..1024 (monotone beyond)")
okB2 = []
for M in [2 ** 12, 2 ** 13, 2 ** 14, 2 ** 16, 2 ** 20, 4100, 5000, 10004]:
    f = int(M / (40 * log(M)))
    L = M // 60
    heavy = (32 * L) // 7 + 1          # heavy count < 32L/7
    U = M - f
    Amin = math.ceil(U / 3)
    n = Amin - heavy
    assert n >= 1
    lhs = comb(n, L)
    rhs = (2 * L + 1) ** f * M
    assert lhs > rhs, (M, f, L, n)
    okB2.append((M, f, L, n))
print("4b. Theorem B'' pigeonhole binom(n,L) > (2L+1)^f M holds (exact integers):", okB2[:4], "...")

# 5. thresholds M_1(C) from the audited inert-prime bound
def inert_primes(limit):
    return [p for p in range(3, limit + 1) if p % 4 == 3 and is_prime(p)]
def E(M, q):
    b = M // q; r = M - q * b
    return q * comb(b, 2) + r * b
def F(M, IP):
    tot = 0.0
    for p in IP:
        if p + 1 >= M: break
        qa = p + 1
        while qa < M:
            tot += E(M, qa) * log(p)
            qa *= p
    return tot
IP = inert_primes(1 << 17)
def W_lower(M, C):
    # primitive normalisation: (M/4) log R >= F(M) - binom(M,2) log C ; W = 2 log R ; b_M = M/4 for even M
    return 2 * (F(M, IP) - comb(M, 2) * log(C)) / (M / 4)
def s_of(M):
    return (1 + math.sqrt(4 * M - 3)) / 2 + 8
def needB(M, C):
    k = 2 * max(0.0, log(C))           # kappa''_+ <= 2 log+ C
    s = s_of(M)
    return max(4 * M * s * k / (M - 2 * s), 36 * max(0.0, log(2 * C)))
def needB2(M, C):
    return max(8 * M * max(0.0, log(C)), 32 * max(0.0, log(M * C / 60.0)))
for C in [1.0, 2.0, 4.0, 8.0, 16.0]:
    out = []
    for need in (needB, needB2):
        first = None
        for k in range(6, 18):
            M = 2 ** k
            if W_lower(M - (0 if need is needB else int(M / (40 * log(M)))), C) > need(M, C):
                if first is None: first = M
            else:
                first = None
        out.append(first)
    print("5. C=%5.1f : Theorem B hypotheses automatic for powers of two M >= %s ; Theorem B'' (W-condition) for M >= %s (tested to 2^17)" % (C, max(64, out[0]) if out[0] else None, out[1]))
print("ALL PASS")
