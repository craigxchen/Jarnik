#!/usr/bin/env python3
"""T3, independent: the 2026 five-divisor #887 construction n=(x^2-1)(x^2-4)(x^2-16)(x^2-49), z^2=56x^2+65.
F0=(x-1)(x-2)(x-4)(x+7), F1=x^4-35x^2+56+2xz, F2=(x^2-4)(x^2-16), F3=(x^2-1)(x^2-16), F4=(x^2-1)(x^2-4).
Checks: F1*conj(F1)=n; split twin u_a = prod_{A ni a} u_A prod_{A not ni a} v_A, u_5 = prod v_A; residues U_a - V_a;
alpha = log H/(-log dist), H=(max/min) lcm; chain deficits Delta_S for S={0, e_a-e_b+e_c-e_d}; w = log lcm/15;
and (stab Prop 3.5) whether Delta_S >= s_S w ever happens (it must not, or the point would have to satisfy x^m = 1)."""
import itertools, math
import flint

def pell_points(xmax):
    pts = set()
    for z0, x0 in [(11, 1), (11, -1), (17, 2), (17, -2), (31, 4), (31, -4), (53, 7), (53, -7)]:
        z, x = z0, x0
        for _ in range(30):
            if 7 < x <= xmax: pts.add((x, z))
            z, x = 15 * z + 112 * x, 2 * z + 15 * x
            if x > xmax: break
    return sorted(pts)

def factor_int(n):
    return {int(p): int(e) for p, e in flint.fmpz(n).factor()[1]} if False else dict((int(p), int(e)) for p, e in flint.fmpz(n).factor())

def vp(n, p):
    v = 0
    while n % p == 0: n //= p; v += 1
    return v

CUTS = [frozenset(A) for r in range(1, 5) for A in itertools.combinations(range(1, 5), r)]

def analyse(x, z):
    F = [(x-1)*(x-2)*(x-4)*(x+7), x**4 - 35*x**2 + 56 + 2*x*z, (x*x-4)*(x*x-16), (x*x-1)*(x*x-16), (x*x-1)*(x*x-4)]
    Fb1 = x**4 - 35*x**2 + 56 - 2*x*z
    n = (x*x-1)*(x*x-4)*(x*x-16)*(x*x-49)
    assert z*z == 56*x*x + 65 and F[1]*Fb1 == n
    g = 0
    for f in F: g = math.gcd(g, f)
    u = [f // g for f in F]
    # label: points 1..4 = F1..F4, point 5 = F0 (base)
    pts = [u[1], u[2], u[3], u[4], u[0]]
    primes = set()
    for r in (1, 2, 4, 7):
        for s in (x - r, x + r):
            primes |= set(factor_int(s).keys())
    prof = {}
    for p in primes:
        a = [vp(q, p) for q in pts]
        if max(a) > 0: prof[p] = a
    # check the tuple is supported on these primes
    for q in pts:
        r = q
        for p in prof: r //= p ** vp(r, p)
        assert r == 1
    uA = {A: 1 for A in CUTS}; vA = {A: 1 for A in CUTS}
    for p, a in prof.items():
        for t in range(1, max(a) + 1):
            S = frozenset(j + 1 for j in range(5) if a[j] >= t)
            if 5 not in S: uA[S] *= p
            else: vA[frozenset(range(1, 6)) - S] *= p
    # twin representation
    for aidx in range(1, 5):
        ua = math.prod(uA[A] for A in CUTS if aidx in A) * math.prod(vA[A] for A in CUTS if aidx not in A)
        assert ua == pts[aidx - 1], 'twin fails'
    assert math.prod(vA.values()) == pts[4]
    U = [math.prod(uA[A] for A in CUTS if a in A) for a in range(1, 5)]
    V = [math.prod(vA[A] for A in CUTS if a in A) for a in range(1, 5)]
    res = [U[i] - V[i] for i in range(4)]
    for i in range(4):
        assert pts[i] * V[i] == pts[4] * U[i]
    lcm = 1
    for q in pts: lcm = lcm * q // math.gcd(lcm, q)
    nA = {A: uA[A] * vA[A] for A in CUTS}
    assert math.prod(nA.values()) == lcm
    H = flint.arb(max(pts)) / min(pts) * lcm
    dist = flint.arb(max(pts) - min(pts)) / min(pts)
    alpha = H.log() / (-dist.log())
    w = math.log(lcm) / 15
    bal = max(abs(math.log(nA[A]) - w) for A in CUTS) / w
    # chain deficits
    best = (0, None); viol = []
    for a, b, c, d in itertools.permutations(range(5), 4):
        if not (a < c and b < d): continue
        m = [0] * 5; m[a] += 1; m[b] -= 1; m[c] += 1; m[d] -= 1
        widthcut = lambda A: abs(sum(m[j - 1] for j in A))
        Wsum = sum(widthcut(A) for A in CUTS)
        s = Wsum - 8  # t = 1 for x^m - 1
        Delta = sum(widthcut(A) * math.log(nA[A]) for A in CUTS) - sum(abs(sum(mi * ai for mi, ai in zip(m, ap))) * math.log(p) for p, ap in prof.items())
        if Delta > best[0]: best = (Delta, (a + 1, b + 1, c + 1, d + 1), s)
        if Delta >= s * w: viol.append((a, b, c, d))
    nonsimple = sum(1 for a in prof.values() if len(set(a)) > 2)
    return dict(x=x, res=res, alpha=float(alpha.mid()), w=w, bal=bal, maxDelta=best, viol=viol, nonsimple=nonsimple, nprimes=len(prof))

if __name__ == '__main__':
    pts = pell_points(4 * 10**10)
    print('Pell points with 7 < x <= 4e10:', len(pts))
    for x, z in pts:
        r = analyse(x, z)
        print(f"x={r['x']:>14d} w={r['w']:6.3f} balance_dev={r['bal']:.3f} alpha={r['alpha']:.4f} residues U_a-V_a={r['res']} "
              f"nonsimple primes={r['nonsimple']}/{r['nprimes']} max chain deficit={r['maxDelta'][0]:.3f} at m={r['maxDelta'][1]} (s={r['maxDelta'][2]}) "
              f"#(Delta>=s*w)={len(r['viol'])}")
