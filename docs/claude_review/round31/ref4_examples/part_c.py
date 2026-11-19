"""(c) Structured families: line arrangements (dual Hesse, Hesse, Moebius-Kantor,
B3, A3), zero-sum parabolas, torsion translates, finite subgroups of PGL_2,
Lagrange interpolation with prescribed cuts.  For each: L, s_h, delta, E and
the slack of (IV.2), plus the Chern side (direct = ledger checked)."""
import itertools, random, sys
import sympy as sp
from front import t, contact_data, iv2_data, greedy_marks, key, report
from engine import chern_side
from collections import Counter

w = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2      # primitive cube root of unity
EXT = sp.sqrt(-3)


def simp(e):
    return sp.nsimplify(sp.expand(e), [sp.sqrt(3)]) if e.has(sp.I) else sp.expand(e)


def lines_to_columns(lines, O, Yf=(3, -7, 2), seed=1):
    """lines: covectors (a,b,c) in P^2.  O: point (vertex of the pencil = fibres).
    New coordinates X' = l1, Z' = l3 vanishing at O, Y' = Yf."""
    O = sp.Matrix(O)
    rnd = random.Random(seed)
    while True:
        l1 = sp.Matrix([rnd.randint(-9, 9) for _ in range(3)])
        l3 = sp.Matrix([rnd.randint(-9, 9) for _ in range(3)])
        # project onto the forms vanishing at O
        l1 = l1 - (l1.dot(O) / O.dot(O)) * O
        l3 = l3 - (l3.dot(O) / O.dot(O)) * O
        Y = sp.Matrix(Yf)
        M = sp.Matrix.vstack(l1.T, Y.T, l3.T)    # new coords = M * old
        if M.det() != 0 and Y.dot(O) != 0:
            break
    Minv = M.inv()
    cols = []
    for l in lines:
        lv = sp.Matrix(l).T * Minv          # l in new coordinates: alpha X' + beta Y' + gamma Z'
        al, be, ga = [simp(z) for z in lv]
        assert be != 0, "line through O"
        cols.append((simp(-(al * t + ga)), be))
    return cols


def run_cols(label, cols, ext=None, N=None, a=None):
    cd = contact_data(cols, ext)
    if N is None:
        N, a = greedy_marks(cd)
    ar = iv2_data(cd, N, a)
    ch = chern_side(cd['m'], cd['finite'], cd['vinf'], cd['degs'])
    lay = Counter()
    for k, (p, e) in N.items():
        for d in range(1, e + 1):
            lay[sum(1 for i in range(cd['m']) if a[i].get(k, 0) >= d)] += p.degree()
    print(f"[{label}] m={cd['m']} layer sizes {dict(lay)}")
    report(label, ch, ar)
    return cd, ch, ar


which = sys.argv[1:] or ['dualhesse', 'hesse', 'mk', 'b3', 'a3', 'parab', 'torsion', 'dihedral', 'octa', 'lagrange']

if 'dualhesse' in which:
    L9 = []
    for k in range(3):
        L9 += [(1, -w ** k, 0), (0, 1, -w ** k), (-w ** k, 0, 1)]
    cols = lines_to_columns(L9, (2, 5, 11))
    run_cols("dual Hesse 9 lines, 12 triple points", cols, EXT)

if 'hesse' in which or 'mk' in which:
    P9 = []
    for k in range(3):
        P9 += [sp.Matrix([0, 1, -w ** k]), sp.Matrix([-w ** k, 0, 1]), sp.Matrix([1, -w ** k, 0])]
    lines = []
    for i, j in itertools.combinations(range(9), 2):
        l = P9[i].cross(P9[j]).applyfunc(simp)
        on = [k for k in range(9) if simp(l.dot(P9[k])) == 0]
        if len(on) == 3 and min(on) == i and sorted(on)[1] == j:
            lines.append((tuple(l), tuple(on)))
    assert len(lines) == 12
    if 'hesse' in which:
        cols = lines_to_columns([l for l, _ in lines], (2, 5, 11))
        run_cols("Hesse 12 lines (greedy layers)", cols, EXT)
    if 'mk' in which:
        sub = [l for l, on in lines if 0 not in on]
        cols = lines_to_columns(sub, (2, 5, 11))
        run_cols("Hesse minus the 4 lines through a point (Moebius-Kantor 8_3)", cols, EXT)

if 'b3' in which:
    B3 = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]
    run_cols("B3 arrangement (9 lines)", lines_to_columns(B3, (2, 5, 11)))
if 'a3' in which:
    A3 = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, -1, 0), (1, 0, -1), (0, 1, -1)]
    run_cols("A3 arrangement (6 lines)", lines_to_columns(A3, (2, 5, 11)))

if 'parab' in which:
    R = [1, 2, 3, 4, -1, -2, -3, -4]
    sets = [S for S in itertools.combinations(R, 4) if sum(S) == 0]
    cols = [(sp.expand(t ** 4 - sp.prod([t - x for x in S])), 1) for S in sets]
    cd = contact_data(cols)
    N, a = {}, [dict() for _ in sets]
    for x in R:
        p = sp.Poly(t - x, t, domain='QQ')
        N[key(p)] = (p, 1)
        for i, S in enumerate(sets):
            if x in S:
                a[i][key(p)] = 1
    run_cols("zero-sum parabolas y = t^4 - prod_S (t-r), R = +-1..+-4", cols, None, N, a)

if 'torsion' in which:
    z6 = -w ** 2
    cols = [(simp(z6 ** i * t ** 2 + t), 1) for i in range(6)]
    run_cols("torsion translates y = zeta_6^i t^2 + t", cols, EXT)
    cols = [(simp((t - w ** i) ** 3), 1) for i in range(3)] + [(simp((t + w ** i) ** 3), 1) for i in range(3)]
    run_cols("Fermat-type y = (t -+ w^i)^3", cols, EXT)

if 'dihedral' in which:
    z6 = -w ** 2
    mats = [((z6 ** k, 0), (0, 1)) for k in range(6)] + [((0, z6 ** k), (1, 0)) for k in range(6)]
    cols = [(simp(a * t + b), simp(c * t + d)) for (a, b), (c, d) in mats]
    run_cols("graphs of the dihedral group D6 < PGL2", cols, EXT)

if 'octa' in which:
    I = sp.I
    gens = [sp.Matrix([[I, 0], [0, 1]]), sp.Matrix([[1, 1], [-1, 1]])]

    def normal(M):
        M = M.applyfunc(sp.expand)
        for z in M:
            if z != 0:
                return (M / z).applyfunc(lambda u: sp.nsimplify(sp.expand(u)))

    G = {tuple(normal(sp.eye(2)))}
    frontier = [sp.eye(2)]
    while frontier:
        nf = []
        for M in frontier:
            for g in gens:
                Nn = normal(g * M)
                k = tuple(Nn)
                if k not in G:
                    G.add(k)
                    nf.append(Nn)
        frontier = nf
    assert len(G) == 24, len(G)
    cols = [(sp.expand(a * t + b), sp.expand(c * t + d)) for a, b, c, d in G]
    run_cols("graphs of the octahedral group S4 < PGL2 (24 sections)", cols, I)

if 'lagrange' in which:
    rnd = random.Random(3)
    Lp = 12
    rs = list(range(1, Lp + 1))
    vals = {r: rnd.randint(-5, 5) for r in rs}
    for m, size in ((10, 6), (16, 6), (16, 8)):
        sets = [tuple(sorted(rnd.sample(rs, size))) for _ in range(m)]
        cols = []
        for S in sets:
            f = sp.interpolate([(r, vals[r]) for r in S], t)
            cols.append((sp.expand(f), 1))
        cd = contact_data(cols)
        N, a = {}, [dict() for _ in sets]
        for r in rs:
            p = sp.Poly(t - r, t, domain='QQ')
            N[key(p)] = (p, 1)
            for i, S in enumerate(sets):
                if r in S:
                    a[i][key(p)] = 1
        run_cols(f"Lagrange: L={Lp}, m={m}, |cut|={size}, f_i = interpolant on its cut", cols, None, N, a)
