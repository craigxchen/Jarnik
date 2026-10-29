#!/usr/bin/env python3
"""T2, independent: gauss Thm 4.1 instance X=-3, m=3 (k=4).  Roots r_T, T subset {1,2,3}:
r_12=-6-4i, r_13=-6-3i, r_23=-6+5i, r_1=-7+6i, r_2=-3-2i, r_3=-4-3i, r_123=i ; H_T(s) = s - r_T.
z_0 = prod conj(H_T), z_x = prod_{T ni x} H_T prod_{T not ni x} conj(H_T), P_x = prod_{T ni x}(s - r_T).
Checks: Im P_x(s) constant; reduced tuples: gcd norm, residues |c_ab|/2, balance, chain deficits for
S={0, e_a-e_b+e_c-e_d} (k=4: s_S = sum_cuts width - 4t = 2), and that Prop 3.5's k=4 analogue never fires."""
import itertools, math, random
import sympy as sp
import flint
sys_path = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/kill_stab_checks'
import sys; sys.path.insert(0, sys_path)
from cox import gm, gc, gn, gsub, ggcd, gdivexact, gval

R = {(1, 2): (-6, -4), (1, 3): (-6, -3), (2, 3): (-6, 5), (1,): (-7, 6), (2,): (-3, -2), (3,): (-4, -3), (1, 2, 3): (0, 1)}
BL = list(R.keys())

def symbolic():
    s = sp.symbols('s', real=True)
    out = []
    for x in (1, 2, 3):
        P = sp.expand(sp.prod([s - (R[T][0] + sp.I * R[T][1]) for T in BL if x in T]))
        out.append(sp.simplify(sp.im(P)))
    return out

def tuple_at(s):
    H = {T: (s - R[T][0], -R[T][1]) for T in BL}
    z0 = (1, 0)
    for T in BL: z0 = gm(z0, gc(H[T]))
    Z = [z0]
    for x in (1, 2, 3):
        z = (1, 0)
        for T in BL: z = gm(z, H[T] if x in T else gc(H[T]))
        Z.append(z)
    return Z, H

def analyse(s):
    Z, H = tuple_at(s)
    g = (0, 0)
    for z in Z: g = ggcd(g, z)
    Zr = [gdivexact(z, g) for z in Z]
    N = gn(Zr[0])
    # residues
    res = []
    for a, b in itertools.combinations(range(4), 2):
        c = gdivexact(gsub(Zr[a], Zr[b]), ggcd(Zr[a], Zr[b]))
        res.append(gn(c))  # |c_ab|^2
    # primes: factor norms of H_T(s)
    primes = set()
    for T in BL:
        for p, e in flint.fmpz(gn(H[T])).factor():
            primes.add(int(p))
    prof = {}
    for p in primes:
        if N % p: continue
        if p % 4 != 1: prof[p] = None; continue
        g0 = 2
        while pow(g0, (p - 1) // 2, p) != p - 1: g0 += 1
        a0 = pow(g0, (p - 1) // 4, p)
        pi = ggcd((p, 0), (a0, 1))
        assert gn(pi) == p
        a = [gval(z, pi) for z in Zr]
        e = gval(Zr[0], pi) + gval(Zr[0], gc(pi))
        prof[p] = (a, e)
    assert all(v is not None for v in prof.values()), 'non-split prime in reduced N'
    # cuts of [4] = {0,1,2,3} named by side A subset {1,2,3} (point 0 = base); layers
    cuts = [frozenset(A) for r in range(1, 4) for A in itertools.combinations((1, 2, 3), r)]
    nA = {A: 1 for A in cuts}
    for p, (a, e) in prof.items():
        for t in range(1, e + 1):
            S = frozenset(j for j in range(4) if a[j] >= t)
            A = S if 0 not in S else frozenset(range(4)) - S
            nA[A] *= p
    assert math.prod(nA.values()) == N
    w = math.log(N) / 7
    bal = max(abs(math.log(nA[A]) - w) for A in cuts) / w
    maxD = 0; fire = 0
    for a, b, c, d in itertools.permutations(range(4), 4):
        m = [0] * 4; m[a] += 1; m[b] -= 1; m[c] += 1; m[d] -= 1
        width = lambda A: abs(sum(m[j] for j in A))
        sS = sum(width(A) for A in cuts) - 4
        D = sum(width(A) * math.log(nA[A]) for A in cuts) - sum(abs(sum(mi * ai for mi, ai in zip(m, ap[0]))) * math.log(p) for p, ap in prof.items())
        maxD = max(maxD, D)
        if D >= sS * w: fire += 1
    nonsimple = sum(1 for ap in prof.values() if len(set(ap[0])) > 2)
    return dict(s=s, gcdnorm=gn(g), maxc2=max(res), w=w, bal=bal, maxD=maxD, fire=fire, nonsimple=nonsimple, np=len(prof))

if __name__ == '__main__':
    print('Im P_x(s) (sympy, s real):', symbolic())
    rnd = random.Random(4)
    S = [10**3 + 7, 10**4 + 1, 123457, 10**6 + 3, 9876543, 10**8 + 7, 987654321, 10**9 + 9] + [rnd.randint(10**5, 10**9) for _ in range(30)]
    worst = dict(gcdnorm=0, maxc2=0, maxD=0, fire=0, bal=0)
    for s in S:
        r = analyse(s)
        for k in ('gcdnorm', 'maxc2', 'maxD', 'fire'):
            worst[k] = max(worst[k], r[k])
        if s > 10**8: worst['bal'] = max(worst['bal'], r['bal'])
        print(f"s={s:>11d} gcd norm={r['gcdnorm']:>9d} max|c_ab|^2={r['maxc2']:>6d} w={r['w']:6.2f} balance={r['bal']:.3f} "
              f"nonsimple={r['nonsimple']}/{r['np']} max chain deficit={r['maxD']:.2f} fires={r['fire']}")
    print('WORST over', len(S), 'values of s:', worst)
