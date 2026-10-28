"""Exact data of the (L_6)-refuting curve E_t (threshold/L8.md Thm B) in the space-quartic model
    E_t ~= { u^2 = kb (a^2 - s h^2),  w^2 = kc (a^2 - tp h^2) }  in P^3[a:u:w:h],
with  b = (-A1b(T) + a u)/(2 A2b(T)),  c = (-A1c(T) + w)/(2 A2c(T)),  T = a^2.
Returns the 25 contact points B_A (A = cluster side, <= 1 anchor) with cluster value v_A in {x1..x4},
and three auxiliary rational points r_j = (a,b,c) = (-x_j, b'_j, c'_j) (j = 1,2,3; not contacts).
Exact arithmetic (Fraction)."""
from fractions import Fraction as Fr
from fam import family, disc_T, check, ev

def coeffs(PQ, T):
    P, Q = PQ
    return [P[k] - T * Q[k] for k in range(3)]  # A0, A1, A2

def other_root(PQ, T, z0):
    A = coeffs(PQ, T)
    return -A[1] / A[2] - z0

def etdata(xs):
    F, D, ok = check(xs, verbose=False)
    assert ok
    xs = F['xs']
    Db, Dc = D['b'], D['c']
    assert Db[0] == 0 and Db[2] != 0 and Dc[2] == 0 and Dc[1] != 0
    kb = Db[2]; s = -Db[1] / Db[2]
    kc = Dc[1]; tp = -Dc[0] / Dc[1]
    mb, mc = F['maps']['b'], F['maps']['c']

    def to_auw(a, b, c):
        T = a * a
        Ab = coeffs(mb, T); Ac = coeffs(mc, T)
        u = (2 * Ab[2] * b + Ab[1]) / a
        w = 2 * Ac[2] * c + Ac[1]
        assert u * u == kb * (a * a - s), 'u'
        assert w * w == kc * (a * a - tp), 'w'
        return (a, u, w)

    def to_abc(a, u, w):
        T = a * a
        Ab = coeffs(mb, T); Ac = coeffs(mc, T)
        return (a, (-Ab[1] + a * u) / (2 * Ab[2]), (-Ac[1] + w) / (2 * Ac[2]))

    pts = []   # (label, side A (frozenset of names), value index j (1..4), (a,u,w))
    aux = []
    for j, x in enumerate(xs, start=1):
        T = x * x
        bp = other_root(mb, T, x); cp = other_root(mc, T, x)
        for sa in (x, -x):
            for bb in (x, bp):
                for cc in (x, cp):
                    S = frozenset(n for n, val in (('a', sa), ('b', bb), ('c', cc)) if val == x)
                    P = to_auw(sa, bb, cc)
                    if j <= 3:
                        if S:
                            pts.append((('A%d' % j,) + tuple(sorted(S)), frozenset({'A%d' % j}) | S, j, P))
                        else:
                            aux.append(P)
                    else:
                        if len(S) >= 2:
                            pts.append((tuple(sorted(S)), S, 4, P))
    assert len(pts) == 25 and len(aux) == 3
    allp = [P for *_, P in pts] + aux
    assert len(set(allp)) == 28, 'points not distinct'
    return dict(xs=xs, kb=kb, s=s, kc=kc, tp=tp, pts=pts, aux=aux, to_abc=to_abc, to_auw=to_auw, F=F)

if __name__ == '__main__':
    import sys
    xs = [Fr(v) for v in sys.argv[1:5]] if len(sys.argv) >= 5 else [Fr(2), Fr(3), Fr(-5), Fr(1)]
    E = etdata(xs)
    print('kb =', E['kb'], ' s =', E['s'], ' kc =', E['kc'], " t' =", E['tp'])
    for lab, A, j, P in E['pts']:
        print(lab, 'v=x%d' % j, [str(c) for c in P])
    print('aux', [[str(c) for c in P] for P in E['aux']])
