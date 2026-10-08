"""Numerical exploration: linear-block polynomial families of five points
on arcs O(R^(2/5)).

Unknowns: gamma_T in C for the 10 edges T of K_5 (block C_T = alpha_T (n+gamma_T)).
Equations (k=1..5, i=1..4):  Im S_k(i) - Im S_k(5) = 0,
   S_k(i) = sum_{T contains i} gamma_T^k.
Normalisation (real translation and scaling of n): sum Re gamma = 0, sum |gamma|^2 = 10.
A solution is degenerate for the pair (i,j) when the multisets
{gamma_il : l != i,j} and {gamma_jl : l != i,j} coincide.
Levenberg-Marquardt from random starts; reports the best residuals and the
degeneracy status of near-solutions.
"""
import sys
from itertools import combinations

import numpy as np

E = list(combinations(range(5), 2))


def residual(v):
    g = v[:10] + 1j * v[10:]
    out = []
    for k in range(1, 6):
        gk = g ** k
        S = [sum(gk[t] for t, T in enumerate(E) if i in T) for i in range(5)]
        for i in range(4):
            out.append((S[i] - S[4]).imag)
    out.append(g.real.sum())
    out.append((abs(g) ** 2).sum() - 10.0)
    return np.array(out)


def jac(v, h=1e-7):
    f0 = residual(v)
    J = np.zeros((len(f0), len(v)))
    for a in range(len(v)):
        w = v.copy()
        w[a] += h
        J[:, a] = (residual(w) - f0) / h
    return f0, J


def lm(v, iters=400):
    lam = 1e-3
    f, J = jac(v)
    cost = f @ f
    for _ in range(iters):
        A = J.T @ J
        g = J.T @ f
        step = np.linalg.solve(A + lam * np.diag(np.diag(A) + 1e-12), -g)
        w = v + step
        fw = residual(w)
        cw = fw @ fw
        if cw < cost:
            v, cost = w, cw
            f, J = jac(v)
            lam = max(lam / 3, 1e-12)
            if cost < 1e-26:
                break
        else:
            lam *= 4
            if lam > 1e12:
                break
    return v, cost


def degeneracy(v, tol=1e-6):
    g = v[:10] + 1j * v[10:]
    deg = []
    for i, j in combinations(range(5), 2):
        rest = [l for l in range(5) if l not in (i, j)]
        a = sorted([g[E.index(tuple(sorted((i, l))))] for l in rest], key=lambda z: (round(z.real, 5), round(z.imag, 5)))
        b = sorted([g[E.index(tuple(sorted((j, l))))] for l in rest], key=lambda z: (round(z.real, 5), round(z.imag, 5)))
        deg.append(all(abs(x - y) < tol for x, y in zip(a, b)))
    minim = min(abs(g.imag))
    return deg, minim


def main(seed=0, starts=300):
    rng = np.random.default_rng(seed)
    results = []
    for s in range(starts):
        v = rng.normal(size=20)
        v, cost = lm(v)
        deg, minim = degeneracy(v)
        results.append((cost, sum(deg), minim, v))
    results.sort(key=lambda t: t[0])
    zero = [r for r in results if r[0] < 1e-20]
    print('starts=%d converged to zero residual: %d' % (starts, len(zero)))
    from collections import Counter
    print('degenerate-pair counts among exact solutions:', Counter(r[1] for r in zero))
    nd = [r for r in zero if r[1] == 0 and r[2] > 1e-4]
    print('fully nondegenerate exact solutions with Im gamma != 0:', len(nd))
    print('smallest residuals among non-converged:', [float('%.3g' % r[0]) for r in results if r[0] >= 1e-20][:8])
    for r in nd[:3]:
        g = r[3][:10] + 1j * r[3][10:]
        print(np.round(g, 6))


if __name__ == '__main__':
    main(*(int(a) for a in sys.argv[1:]))


def inspect(seed=0, starts=60):
    rng = np.random.default_rng(seed)
    from collections import Counter
    c = Counter()
    ex = {}
    for s in range(starts):
        v = rng.normal(size=20)
        v, cost = lm(v)
        if cost > 1e-20:
            continue
        g = v[:10] + 1j * v[10:]
        nreal = int(sum(abs(g.imag) < 1e-6))
        deg, _ = degeneracy(v)
        key = (nreal, sum(deg))
        c[key] += 1
        ex.setdefault(key, g)
    print('(number of real gammas, degenerate pairs):', c)
    for key, g in ex.items():
        print(key, [complex(round(z.real, 4), round(z.imag, 4)) for z in g])


def generic(g, tol=1e-4):
    if min(abs(g.imag)) < tol:
        return False
    for a in range(10):
        for b in range(10):
            if a != b and (abs(g[a] - g[b]) < tol or abs(g[a] - g[b].conjugate()) < tol):
                return False
    return True


def survey(seed=0, starts=200):
    rng = np.random.default_rng(seed)
    exact = gen = 0
    best_generic_cost = None
    for s in range(starts):
        v = rng.normal(size=20)
        v, cost = lm(v)
        g = v[:10] + 1j * v[10:]
        if cost < 1e-20:
            exact += 1
            if generic(g):
                gen += 1
                print('GENERIC EXACT SOLUTION', [complex(round(z.real, 6), round(z.imag, 6)) for z in g])
    print('seed=%d starts=%d exact=%d generic_exact=%d' % (seed, starts, exact, gen))
