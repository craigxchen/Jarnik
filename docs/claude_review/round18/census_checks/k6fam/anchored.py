"""(1) Open conditions of U at t0; (2) anchored (L_6) form of actual configurations on E_t0 via Lemma 5.2:
P_0=(1,0) for the anchor x1 (sent to infinity), P_j ~ (D X_j, Y_j), D=lcm|Y_j|; report max|det(P0,Pj)| and the
residues t_ij = det(P_i,P_j)/prod_{T contains i,j} n_T with n_T the gcd-extracted contact integers."""
from fractions import Fraction as Fr
from math import gcd, log
from itertools import combinations
from fam import family, fder, disc_T, other_root, ev
from ell import quotient_curve, nontorsion
from lift import lift_T
from profile import profile, sides, mobius_to_std, row, det

def U_conditions(xs):
    F = family(xs); maps = F['maps']; e1, e2, e3, e4 = F['e']
    Db, Dc = disc_T(maps['b']), disc_T(maps['c'])
    s = -Db[1] / Db[2]; t = -Dc[0] / Dc[1]
    conds = {
        'x distinct, x_j != +-x_k, nonzero': all(x != 0 for x in xs) and all(x != y and x != -y for x, y in combinations(xs, 2)),
        'e1 != 0': e1 != 0, 'e3^2-4e2e4 != 0': e3 * e3 - 4 * e2 * e4 != 0,
        'lam_b != lam_c (and != 0)': F['lam']['b'] not in (0, F['lam']['c']) and F['lam']['c'] != 0,
        'Db = T(T-s) shape, s not in {0,inf}': Db[0] == 0 and Db[2] != 0 and s != 0,
        'Dc linear, t not in {0,inf}': Dc[2] == 0 and Dc[1] != 0 and t != 0,
        's != t': s != t,
        'x_j^2 not in {s,t}': all(x * x not in (s, t) for x in xs),
        'Q_lam(x_j) != 0': all(ev(maps[q][1], x) != 0 for q in 'abc' for x in xs),
        'derivatives distinct & nonzero at x_j': all(len({fder(maps[q], x) for q in 'abc'}) == 3 and all(fder(maps[q], x) != 0 for q in 'abc') for x in xs),
    }
    return conds, s, t

def anchored_form(xs, pt):
    d, b, Dm = profile(xs, pt)          # contact integers for sides C (<=1 anchor), anchors at inf,0,1
    x1, x2, x3 = xs[:3]
    M = mobius_to_std(x1, x2, x3)
    pts = [None, Fr(0), Fr(1)] + [M(z) for z in pt]
    V = [row(p) for p in pts]          # V[0] = (1,0) is the anchor x1
    Ys = [abs(V[j][1]) for j in range(1, 6)]
    Dl = 1
    for y in Ys:
        Dl = Dl * y // gcd(Dl, y)
    P = [(1, 0)]
    for j in range(1, 6):
        X, Y = V[j][0] * Dl, V[j][1]
        g = gcd(X, Y); P.append((X // g, Y // g))
    # blocks n_T for T subset {1..5} (side not containing anchor 0): T = C if 0 not in C else complement
    SC = sides(); n = {}
    for C in SC:
        T = C if 0 not in C else frozenset(range(6)) - C
        n[T] = n.get(T, 1) * d[C]
    Y = [abs(det(P[0], P[j])) for j in range(1, 6)]
    tmax = 0
    for i, j in combinations(range(1, 6), 2):
        pr = 1
        for T, v in n.items():
            if i in T and j in T:
                pr *= v
        dd = det(P[i], P[j])
        assert dd % pr == 0
        tmax = max(tmax, abs(dd // pr))
    return max(Y), tmax, n

if __name__ == '__main__':
    xs = [Fr(2), Fr(3), Fr(-5), Fr(1)]
    conds, s, t = U_conditions(xs)
    print('t0 =', [str(x) for x in xs], ' s =', s, ' t =', t)
    for kk, v in conds.items():
        print('   %-45s %s' % (kk, v))
    F, C, g, pts = quotient_curve(xs)
    k = g[3]
    # non-torsion of 2*Q1 (image of P - O with O=(x1,x1,x1), P=(x1,b1',x1))
    print('Q1 non-torsion:', nontorsion(C, pts[0]), ' 2Q1 non-torsion:', nontorsion(C, C.add(pts[0], pts[0])))
    Q1, Q2 = pts[0], pts[1]
    Dl = C.add(Q2, (Q1[0], -Q1[1])); R = Q1
    for N in range(1, 15):
        R = C.add(R, Dl)
        L = lift_T(F, R[0] / k)
        try:
            Ymax, tmax, n = anchored_form(xs, L[0])
        except ZeroDivisionError:
            continue
        logs = [log(v) for T, v in n.items()]
        w = sum(logs) / len(logs)
        print('N=%2d  anchored: max|Y_j| = %d, max|t_ij| = %d, 25 blocks: min/w=%.3f max/w=%.3f (w=%.1f)' % (N, Ymax, tmax, min(logs) / w, max(logs) / w, w))
