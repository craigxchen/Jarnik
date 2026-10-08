"""Exact checks of the general signed-character defect identity (Theorem A of walsh-extraction.md).

V_c = W/2 - kappa/2 - (1/(2h)) sum_{x!=y in P} c_x c_y u_xy + Delta_c,
Delta_c = (1/2) sum_j w_j |T_j| (h-|T_j|)/h,   u_xy = kappa - G_xy,
for every balanced +-1 character c on h rows, and its general-amplitude version with Lambda >= max|T_j|.
Part 2 realizes the height on literal Gaussian integers: prod z^c = unit * beta/conj(beta), Norm(beta)=exp(V_c).
"""
import random, math, sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gauss_common import *

random.seed(20261008)

def identity_check(S, w, kappa, c):
    M = len(S); K = len(w)
    W = sum(w)
    T = [sum(c[x] * S[x][j] for x in range(M)) for j in range(K)]
    V = sum(Fr(w[j], 2) * abs(T[j]) for j in range(K))
    G = [[sum(w[j] * S[x][j] * S[y][j] for j in range(K)) for y in range(M)] for x in range(M)]
    Lam = max(max(abs(t) for t in T), 1)
    c2 = sum(ci * ci for ci in c)
    sccu = sum(c[x] * c[y] * (kappa - G[x][y]) for x in range(M) for y in range(M) if x != y)
    Delta = sum(Fr(w[j] * abs(T[j]) * (Lam - abs(T[j])), 2 * Lam) for j in range(K))
    rhs = Fr(c2, 2 * Lam) * (W - kappa) - Fr(1, 2 * Lam) * sccu + Delta
    assert V == rhs, (V, rhs)
    assert Delta >= 0
    return V, Delta, sccu, Lam

count = 0; pm_count = 0
for trial in range(3000):
    M = random.randint(4, 12); K = random.randint(3, 20)
    S = [[random.choice((1, -1)) for _ in range(K)] for _ in range(M)]
    w = [random.randint(1, 50) for _ in range(K)]
    kappa = Fr(random.randint(-200, 200), random.randint(1, 9))
    if trial % 2 == 0:
        h = 2 * random.randint(1, M // 2)
        P = random.sample(range(M), h)
        signs = [1] * (h // 2) + [-1] * (h // 2); random.shuffle(signs)
        c = [0] * M
        for x, s in zip(P, signs): c[x] = s
        V, Delta, sccu, Lam = identity_check(S, w, kappa, c)
        # +-1 balanced case with Lambda = h : V = W/2 - kappa/2 - sccu/(2h) + Delta_h
        W = sum(w)
        T = [sum(c[x] * S[x][j] for x in range(M)) for j in range(K)]
        Dh = sum(Fr(w[j] * abs(T[j]) * (h - abs(T[j])), 2 * h) for j in range(K))
        assert V == Fr(W, 2) - kappa / 2 - sccu / (2 * h) + Dh
        pm_count += 1
    else:
        c = [random.randint(-3, 3) for _ in range(M)]
        s = sum(c); c[0] -= s
        if all(ci == 0 for ci in c): continue
        identity_check(S, w, kappa, c)
    count += 1
print("Part 1: exact identity on", count, "random weighted profiles (", pm_count, "balanced +-1 cases)")

# four-row specialisation: Delta = (W - T1234)/4 for c=(1,1,-1,-1)
for trial in range(500):
    K = random.randint(3, 15)
    S = [[random.choice((1, -1)) for _ in range(K)] for _ in range(4)]
    w = [random.randint(1, 30) for _ in range(K)]
    c = [1, 1, -1, -1]
    T = [sum(c[x] * S[x][j] for x in range(4)) for j in range(K)]
    Dh = sum(Fr(w[j] * abs(T[j]) * (4 - abs(T[j])), 8) for j in range(K))
    T1234 = sum(w[j] * S[0][j] * S[1][j] * S[2][j] * S[3][j] for j in range(K))
    assert Dh == Fr(sum(w) - T1234, 4)
print("Part 1b: four-row specialisation Delta=(W-T1234)/4 on 500 profiles")

# Part 2: literal Gaussian realisation
def paley_matrix(q):
    chi = {}
    for a in range(q): chi[a] = 0
    for a in range(1, q): chi[(a * a) % q] = 1
    def ch(a):
        a %= q
        if a == 0: return 0
        return 1 if chi[a] == 1 else -1
    M = q + 1
    H = [[0] * M for _ in range(M)]
    for j in range(M): H[0][j] = 1
    for x in range(M): H[x][0] = 1
    for x in range(1, M):
        for j in range(1, M):
            H[x][j] = -1 if x == j else -ch((x - 1) - (j - 1))
    return H

def walsh_matrix(t):
    M = 2 ** t
    return [[(-1) ** bin(x & a).count("1") for a in range(M)] for x in range(M)]

def flipped_profile(H, b):
    M = len(H); labels = list(range(1, M))
    cols = []
    for a in labels:
        for _ in range(b): cols.append([H[x][a] for x in range(M)])
    # one distinct flip per row, capacity two per label
    order = list(range(M)); random.shuffle(order)
    used = set(); slot = 0
    for idx, x in enumerate(order):
        a = labels[idx % len(labels)]
        cpos = (a - 1) * b + (idx // len(labels))
        cols[cpos][x] *= -1
    return [[cols[j][x] for j in range(len(cols))] for x in range(M)]

gchecks = 0
for name, H, b in [("paley12", paley_matrix(11), 5), ("walsh8", walsh_matrix(3), 5), ("paley20", paley_matrix(19), 5)]:
    S = flipped_profile(H, b)
    M = len(S); K = len(S[0])
    primes = split_primes(K + 3)[3:]
    pis = [gauss_prime(p) for p in primes]
    orient = [random.choice((1, -1)) for _ in range(K)]
    g = (2, 1)  # common content 2+i
    units = [random.choice(UNITS) for _ in range(M)]
    def factor(j, s):
        pi = pis[j] if orient[j] * s == 1 else gconj(pis[j])
        return pi
    z = []
    for x in range(M):
        v = gmul(g, units[x])
        for j in range(K): v = gmul(v, factor(j, S[x][j]))
        z.append(v)
    N = gnorm(z[0]); assert all(gnorm(v) == N for v in z)
    for trial in range(60):
        h = random.choice([2, 4, 6, 8])
        P = random.sample(range(M), h)
        signs = [1] * (h // 2) + [-1] * (h // 2); random.shuffle(signs)
        c = [0] * M
        for x, s in zip(P, signs): c[x] = s
        T = [sum(c[x] * S[x][j] for x in range(M)) for j in range(K)]
        assert all(t % 2 == 0 for t in T)
        v = [t // 2 for t in T]
        Pn = (1, 0); Qd = (1, 0)
        for x in range(M):
            if c[x] > 0: Pn = gmul(Pn, gpow(z[x], c[x]))
            if c[x] < 0: Qd = gmul(Qd, gpow(z[x], -c[x]))
        beta = (1, 0)
        for j in range(K):
            if v[j] != 0:
                base = factor(j, 1 if v[j] > 0 else -1)
                beta = gmul(beta, gpow(base, abs(v[j])))
        lhs = gmul(Pn, gconj(beta))
        ok = any(gmul(gmul(u, Qd), beta) == lhs for u in UNITS)
        assert ok
        Vc = sum(abs(v[j]) * math.log(primes[j]) for j in range(K))
        assert abs(math.log(gnorm(beta)) - Vc) < 1e-9 * max(1, Vc)
        # identity with real weights
        w = [math.log(p) for p in primes]; W = sum(w)
        kappa = 0.37
        G = [[sum(w[j] * S[a][j] * S[bb][j] for j in range(K)) for bb in range(M)] for a in range(M)]
        sccu = sum(c[a] * c[bb] * (kappa - G[a][bb]) for a in range(M) for bb in range(M) if a != bb)
        Dh = sum(w[j] * abs(T[j]) * (h - abs(T[j])) / (2 * h) for j in range(K))
        assert abs(Vc - (W / 2 - kappa / 2 - sccu / (2 * h) + Dh)) < 1e-8 * W
        gchecks += 1
print("Part 2: literal Gaussian heights and identity on", gchecks, "characters (Paley12/Walsh8/Paley20 flipped, content 2+i, random units)")
print("ALL PASS")
