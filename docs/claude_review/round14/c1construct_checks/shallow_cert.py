# Theorem S (exact check): coherent labels with root-of-unity parity repair on the shallow classes
# k(q) = v_2(m_q) <= E0, with auxiliary Gaussian primes zeta_j whose norms l_j are = 1 mod 2^(E0+1)
# (so no l_j is a shallow modulus; the design is defined at EVERY shallow q <= X, any X).
#   t_j(Q) = (zeta_j/conj zeta_j) * f_Q^{eps_j(q)},  eps_j(q) = [(l_j/q) != (p_j/q)],  f_Q of order 2^{k(q)}.
# Checks, at every shallow prime q <= X (q not dividing N) and every level Q = q^a <= X:
#   (1) parity: t_j is a non-square in T_Q iff (p_j/q) = -1  [=> realisable with Norm r_j = p_j, sharp Prop 2.1(4)];
#       level-1 realisation found explicitly for q <= QREAL;
#   (2) for every v with 1 <= ||v||_1 <= KMAX and every unit s: collision at q^a  =>  q^a | Im(gamma_v^(2^k(q)));
#   (3) total shallow demand  sum_q a_{q,s}(v) log q  <=  (2^E0 - 2) sum_j |v_j| log l_j.
import sys, math, random, itertools
sys.path.insert(0, '.')
from lib14 import *
from sympy import primerange, factorint
import lib10

def two_sylow_gen(Q, q, k, rng):
    m = mq(q) * (Q // q)
    while True:
        a, b = rng.randrange(Q), rng.randrange(Q)
        if (a * a + b * b) % q == 0: continue
        t = tor((a, b), Q)
        f = mpow(t, m // (2 ** k), Q)
        if mpow(f, 2 ** (k - 1), Q) != (1, 0):
            return f

def is_nonsquare(t, Q, q):
    m = mq(q) * (Q // q)
    return mpow(t, m // 2, Q) != (1, 0)

def run(ps, E0, X, KMAX, QREAL, seed):
    rng = random.Random(seed)
    r = len(ps)
    ls = [l for l in primerange(3, 10**6) if l % (2 ** (E0 + 1)) == 1][:r]
    zs = [gauss_prime(l) for l in ls]
    N_primes = set(ps)
    vs = lib10.monos(r, KMAX, half=True)
    demand = {(v, s): 0.0 for v in vs for s in range(4)}
    nmod = 0; nlev = 0; ncoll = 0
    for q in primerange(3, X + 1):
        if q in N_primes or kq(q) > E0: continue
        nmod += 1
        k = kq(q)
        a_top = 1
        while q ** (a_top + 1) <= X: a_top += 1
        Qtop = q ** a_top
        ftop = two_sylow_gen(Qtop, q, k, rng) if k > 2 else (0, 1)   # k=2: f = i
        ttop = []
        for j in range(r):
            eps = 1 if (legendre(ls[j], q) != legendre(ps[j], q)) else 0
            t = tor(zs[j], Qtop)
            if eps: t = mmul(t, ftop, Qtop)
            ttop.append(t)
        for a in range(1, a_top + 1):
            Q = q ** a; nlev += 1
            f = (ftop[0] % Q, ftop[1] % Q)
            assert mpow(f, 2 ** k, Q) == (1, 0) and (k == 0 or mpow(f, 2 ** (k - 1), Q) != (1, 0))
            ts = []
            for j in range(r):
                t = (ttop[j][0] % Q, ttop[j][1] % Q)   # lower levels are reductions of the top level
                # (1) parity
                assert is_nonsquare(t, Q, q) == (legendre(ps[j], q) == -1), (q, a, j)
                ts.append(t)
            if a == 1 and q <= QREAL:
                for j in range(r):
                    ok = False
                    for x in range(q):
                        for y in range(q):
                            if (x * x + y * y) % q == ps[j] % q and (x * x + y * y) % q != 0:
                                if tor((x, y), q) == ts[j]:
                                    ok = True; break
                        if ok: break
                    assert ok, ("realise", q, j)
            for v in vs:
                tv = (1, 0)
                for j in range(r):
                    if v[j]:
                        tv = mmul(tv, mpow(ts[j] if v[j] > 0 else minv(ts[j], Q), abs(v[j]), Q), Q)
                for s in range(4):
                    if tv == unit_mod(s, Q):
                        ncoll += 1
                        # (2) certificate
                        g = Zv(zs, v)
                        im = gpow(g, 2 ** k)[1]
                        assert im != 0 and im % Q == 0, (q, a, v, s)
                        demand[(v, s)] += math.log(q)      # level a adds log q on top of level a-1
    worst = 0.0; worst_v = None
    for (v, s), d in demand.items():
        bound = (2 ** E0 - 2) * sum(abs(t) * math.log(l) for t, l in zip(v, ls))
        assert d <= bound + 1e-9, (v, s, d, bound)
        rat = d / bound
        if rat > worst: worst, worst_v = rat, (v, s, d, bound)
    print(f"r={r} E0={E0} X={X} KMAX={KMAX}: zeta norms {ls}; shallow moduli {nmod}, levels {nlev}, collisions {ncoll}")
    print(f"   parity OK at every level; level-1 realisation found for all q <= {QREAL}; every collision certified by Im(gamma^(2^k))")
    print(f"   max ratio (shallow demand)/((2^E0-2) sum|v_j| log l_j) = {worst:.3f} at v={worst_v[0]}, s={worst_v[1]} (demand {worst_v[2]:.2f}, bound {worst_v[3]:.2f})")

def legendre(a, q):
    return lib10.legendre_symbol(a % q, q) if a % q else 0

run([10009, 10037, 10069, 10093, 10133, 10141], 3, 3000, 3, 150, 1)
run([10009, 10037, 10069, 10093, 10133, 10141, 10169, 10177], 4, 2000, 2, 100, 2)
