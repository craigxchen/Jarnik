"""Exact verification of explicit integer two-class families (homogeneous polynomial families).
For each family: list of binary forms y_f(alpha, beta); nodes; we evaluate at coprime (alpha, beta),
compute q2 exactly (tc.q2), the diameter, the ratio, and the per-node check that the values occupy
two classes modulo each node form value."""
from math import gcd, log
from tc import q2, split_gcds

def fam_parallelogram(a, b):        # B=4, E=2 : nodes alpha, beta (d=1)
    return [0, a, b, a + b], [a, b]

def fam5(a, b):                     # B=5, E=2 (d=2): nodes a, b, a+b, a+2b
    return [0, a * (a + b), b * (a + 2 * b), a * (a + b) + b * (a + 2 * b), (a + b) * (a + 2 * b)], [a, b, a + b, a + 2 * b]

def fam6(a, b):                     # B=6, E=2 (d=2): nodes a, b, a-b, a-2b ; symmetric f <-> h - f
    return [0, (a - b) * (a - 2 * b), a * (a - 2 * b), -2 * b * (a - b), -a * b, a * a - 4 * a * b + 2 * b * b], [a, b, a - b, a - 2 * b]

def cube(ms):                       # B=2^n, E=n/(n-1): a_j = M/m_j
    M = 1
    for m in ms: M *= m
    a = [M // m for m in ms]
    pts = []
    for mask in range(1 << len(ms)):
        pts.append(sum(a[j] for j in range(len(ms)) if mask >> j & 1))
    return pts, ms

def report(name, pts, nodes):
    pts = sorted(set(pts))
    L = pts[-1] - pts[0]
    q = q2(pts)
    prodn = 1
    for n in nodes: prodn *= abs(n)
    ok = all(len({p % abs(n) for p in pts}) <= 2 for n in nodes if abs(n) > 1)
    print(f"{name}: B={len(pts)} logL={log(L):.2f} log q2={log(q):.2f} ratio={log(q)/log(L):.4f} "
          f"two classes mod each node value: {ok}; q2 divisible by prod|nodes|: {q % prodn == 0}")

if __name__ == "__main__":
    for X in [10**3, 10**6, 10**12, 10**24]:
        a, b = 3 * X + 1, X
        report(f"parallelogram X={X:.0e}", *fam_parallelogram(a, b))
        report(f"fam5 X={X:.0e}", *fam5(a, b))
        report(f"fam6 X={X:.0e}", *fam6(a, b))
    import sympy
    for n in [3, 4, 5]:
        for bits in [20, 60]:
            ms = [sympy.prevprime(2 ** bits - 1000 * j) for j in range(n)]
            report(f"cube n={n} bits={bits}", *cube(ms))
