"""(d) Characteristic-2 trace construction: own GF(2^d) arithmetic.
Checks the incidence/divisibility claims, computes both sides of (IV.2) and the
quantity 3 e_orb - (K+B)^2 of Astra's surface (which Langer's inequality would
force to be >= 0 in characteristic zero)."""
import sys
from engine import chern_side, iv2

POLY = {3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011, 7: 0b10000011}


def run(d):
    q = 1 << d
    mod = POLY[d]

    def mul(x, y):
        r = 0
        while y:
            if y & 1:
                r ^= x
            y >>= 1
            x <<= 1
            if x & q:
                x ^= mod
        return r

    def pw(x, n):
        r = 1
        while n:
            if n & 1:
                r = mul(r, x)
            x = mul(x, x)
            n >>= 1
        return r

    inv = lambda x: pw(x, q - 2)
    tr = lambda x: [0, 1].index(sum_tr(x))

    def sum_tr(x):
        s, y = 0, x
        for _ in range(d):
            s ^= y
            y = mul(y, y)
        assert s in (0, 1)
        return s

    def f_coeffs(a):
        # L_a(t) = a^{-q/2} sum_j (a t)^{2^j};  f_a = (t^{q/2} - L_a)/t ; char 2: minus = plus
        c = [0] * (q // 2 + 1)
        ai = pw(inv(a), q // 2)
        for j in range(d):
            e = 1 << j
            c[e] ^= mul(ai, pw(a, e))
        c[q // 2] ^= 1
        assert c[0] == 0
        f = c[1:]
        while f and f[-1] == 0:
            f.pop()
        return f

    def ev(f, x):
        r = 0
        for co in reversed(f):
            r = mul(r, x) ^ co
        return r

    marks = list(range(1, q))
    m = len(marks)
    F = [f_coeffs(a) for a in marks]
    degs = [len(f) - 1 for f in F]
    pts = list(range(q))  # all of F_q: finite points where roots may lie
    finite = []
    # valuations: check D_ab has all roots simple and in F_q^*
    Dd = {}
    for i in range(m):
        for j in range(i + 1, m):
            n = max(len(F[i]), len(F[j]))
            D = [(F[i][k] if k < len(F[i]) else 0) ^ (F[j][k] if k < len(F[j]) else 0) for k in range(n)]
            while D and D[-1] == 0:
                D.pop()
            Dd[(i, j)] = D
            roots = [x for x in pts if ev(D, x) == 0]
            common = [x for x in range(1, q) if sum_tr(mul(marks[i], x)) == 0 and sum_tr(mul(marks[j], x)) == 0]
            assert roots == common, (i, j)
            assert len(D) - 1 == len(common) == q // 4 - 1
            # roots simple: derivative nonzero at roots (formal derivative in char 2)
            Dp = [D[k] if k % 2 == 1 else 0 for k in range(1, len(D))]
            assert all(ev(Dp, x) != 0 for x in roots)
    for x in pts:
        v = [[0] * m for _ in range(m)]
        any_ = False
        for (i, j), D in Dd.items():
            if ev(D, x) == 0:
                v[i][j] = v[j][i] = 1
                any_ = True
        if any_:
            finite.append((1, v))
    vinf = [[0] * m for _ in range(m)]
    for (i, j), D in Dd.items():
        vinf[i][j] = vinf[j][i] = degs[i] + degs[j] - (len(D) - 1)
    ch = chern_side(m, finite, vinf, degs)
    # IV.2 data: N = prod_{x != 0}(t - x), n_a = prod_{Tr(ax)=0}(t-x)
    L = q - 1
    layers = [(1, sum(1 for a in marks if sum_tr(mul(a, x)) == 0)) for x in range(1, q)]
    alpha = [sum(1 for x in range(1, q) if sum_tr(mul(a, x)) == 0) for a in marks]
    degD = {k: len(D) - 1 for k, D in Dd.items()}
    E = {k: 0 for k in Dd}   # verified above: every common marked x is a root
    ar = iv2(m, L, layers, alpha, degD, E)
    print(f"q={q} m={m} deg f={sorted(set(degs))} S={ch['S']} I={ch['I']} infsteps={ch['steps']} big={sorted(set(ch['big']))}x{len(ch['big'])} "
          f"unblown={ch['unblown']} | 3e_orb-(K+B)^2={ch['Q']} (={float(ch['Q']):.3f}) ledger={ch['ledger']} || "
          f"IV.2: L={ar['L']} delta={ar['delta']} lhs={ar['lhs']} rhs={ar['rhs']} slack={ar['slack']}", flush=True)


for d in map(int, sys.argv[1:]):
    run(d)
