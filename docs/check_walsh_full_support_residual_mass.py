#!/usr/bin/env python3
"""Finite algebra audit; fixtures are not asserted to satisfy arc phases."""
import math
import random


def parity(n):
    return bin(n).count('1') % 2


def subspaces(M):
    spaces = {frozenset([0])}
    frontier = list(spaces)
    while frontier:
        V = frontier.pop()
        for x in range(1, M):
            if x not in V:
                U = V | frozenset(v ^ x for v in V)
                if U not in spaces:
                    spaces.add(U)
                    frontier.append(U)
    return spaces


def close(x, y):
    assert abs(x-y) < 1e-8 * max(1, abs(x), abs(y)), (x, y)


def audit(M, rng):
    b = 5
    labels = [a for a in range(1, M) for _ in range(b)]
    weights = [rng.randint(2, 100) for _ in labels]
    flips = rng.sample(range(len(labels)), M)
    assigned = [labels[j] for j in flips]
    f = [weights[j] for j in flips]
    Wa = [0] * M
    Fa = [0] * M
    for a, w in zip(labels, weights):
        Wa[a] += w
    for a, w in zip(assigned, f):
        Fa[a] += w
    W, F = sum(Wa), sum(f)
    C, D = 1.7, 3.2
    A = lambda h: math.log(8)-2*math.log(h*C)
    q = [Wa[a]-4*Fa[a]/M for a in range(M)]
    low = (W+D-2*F+2*A(M))/M
    r = [0] + [q[a]-low for a in range(1, M)]
    close(sum(r), W/M+(2-6/M)*F-(1-1/M)*(D+2*A(M)))
    for a in range(1, M):
        close(M*Wa[a]/2+F-2*Fa[a], M*q[a]/2+F)
    count = 0
    for V in subspaces(M):
        h = len(V)
        if h < 2:
            continue
        chars = {}
        ordered = sorted(V)
        for a in range(1, M):
            pattern = tuple(parity(a & v) for v in ordered)
            if any(pattern):
                chars.setdefault(pattern, a)
        cosets = []
        unseen = set(range(M))
        while unseen:
            x0 = min(unseen)
            P = {x0 ^ v for v in V}
            unseen -= P
            cosets.append((x0, P))
        for pattern, a0 in chars.items():
            cell = [a for a in range(1, M)
                    if tuple(parity(a & v) for v in ordered) == pattern]
            assert len(cell) == M//h
            Rc, Fc = sum(r[a] for a in cell), sum(Fa[a] for a in cell)
            gaps = []
            for x0, P in cosets:
                fp = sum(f[x] for x in P)
                fm = sum(f[x] for x in P if assigned[x] in cell)
                height = h*sum(Wa[a] for a in cell)/2+fp-2*fm
                residual_gap = h*Rc/2+2*h*Fc/M+fp-2*fm-F-2*math.log(M/h)
                close(residual_gap, height-W/2-D/2-A(h))
                c = [(1 if parity(a0 & (x ^ x0)) == 0 else -1)
                     if x in P else 0 for x in range(M)]
                direct = 0
                for j, (label, weight) in enumerate(zip(labels, weights)):
                    twice = sum(c[x]*(1 if parity(label & x) == 0 else -1)
                                *(-1 if flips[x] == j else 1) for x in range(M))
                    assert twice % 2 == 0
                    direct += abs(twice//2)*weight
                close(height, direct)
                gaps.append(residual_gap)
                count += 1
            close(sum(gaps)/len(gaps), h*Rc/2+h*F/M-F-2*math.log(M/h))
    return count


if __name__ == '__main__':
    rng = random.Random(20260915)
    count = sum(audit(M, rng) for M in [8, 8, 16])
    print('PASS: {} exact coset heights and residual-gap/average identities'.format(count))
