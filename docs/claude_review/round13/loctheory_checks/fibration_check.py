"""Check of the string (fibration) theorem on small column-cut ternary boxes.

Box: A = {k in {-1,0,1}^d : |sum_{x in S_p} k_x| <= L_p for all p}, cuts S_p.
Atoms X = Venn cells of the S_p (including the outside cell).  For each atom of size u the Lefschetz strings of
k[x]/(x^3) under sigma_X are: for i = 0..u, (h_i - h_{i-1}) strings of radius r = u - i starting in degree i.

Predictions checked (all exact over Q):
 (P1) H_{R/J(A)} computed directly from A (reduced monomials, exponents <= 2)
 (P2) sum_s t^{b_s} H_{Bx_s cap P}            (conjectural 'fibred Hilbert function')
 (P3) sum_s t^{b_s} H_{S_a/(J(Bx_s)+J(P))}     (= Hilbert function of R/((x^3)+pi^*J(P)), theorem)
and BIP(s): J(Bx_s cap P) = J(Bx_s) + J(P) for all string boxes.
"""
import itertools, sys
from math import comb
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from lfi import hilbert, monos, leading_forms, ideal_dim_in_degree

def trinom(u):
    h = [1]
    for _ in range(u):
        nh = [0] * (len(h) + 2)
        for a, c in enumerate(h):
            nh[a] += c; nh[a + 1] += c; nh[a + 2] += c
        h = nh
    return h

def strings(u):
    h = trinom(u)
    out = []
    for i in range(u + 1):
        m = h[i] - (h[i - 1] if i > 0 else 0)
        if m:
            out.append((i, u - i, m))  # (shift, radius, multiplicity)
    return out

def box_points(d, cuts, Ls):
    pts = []
    for k in itertools.product((-1, 0, 1), repeat=d):
        if all(abs(sum(k[x] for x in S)) <= L for S, L in zip(cuts, Ls)):
            pts.append(k)
    return pts

def atoms_of(d, cuts):
    pat = {}
    for x in range(d):
        key = tuple(int(x in S) for S in cuts)
        pat.setdefault(key, []).append(x)
    return pat  # pattern -> list of coords

def atom_set(atoms, radii, cuts_pat, Ls):
    """lattice points t in prod [-r_X, r_X] with |sum_{X: pat_X[p]=1} t_X| <= L_p"""
    keys = list(atoms)
    pts = []
    for t in itertools.product(*[range(-radii[X], radii[X] + 1) for X in keys]):
        ok = True
        for p, L in enumerate(Ls):
            s = sum(tv for tv, X in zip(t, keys) if X[p] == 1)
            if abs(s) > L:
                ok = False; break
        if ok:
            pts.append(t)
    return pts

def addH(tot, H, shift):
    for i, v in enumerate(H):
        while len(tot) <= i + shift:
            tot.append(0)
        tot[i + shift] += v

def run(d, cuts, Ls, check_bip=True):
    A = box_points(d, cuts, Ls)
    HA = hilbert(A, 2 * d, cap=2)
    atoms = atoms_of(d, cuts)
    keys = list(atoms)
    D = len(keys)
    full = {X: len(atoms[X]) for X in keys}
    P = atom_set(atoms, full, None, Ls)
    # J(P) generators: leading forms in degrees up to rho(P)+1
    HP = hilbert(P, 10 ** 6)
    rP = len(HP) - 1
    JP = []
    for i in range(rP + 2):
        top, basis = leading_forms(P, i)
        for v in basis:
            JP.append({e: c for e, c in zip(top, v) if c})
    pred2, pred3 = [], []
    bip_fail = []
    for combo in itertools.product(*[strings(full[X]) for X in keys]):
        shift = sum(c[0] for c in combo)
        radii = {X: c[1] for X, c in zip(keys, combo)}
        mult = 1
        for c in combo:
            mult *= c[2]
        Q = atom_set(atoms, radii, None, Ls)
        HQ = hilbert(Q, 10 ** 6)
        for _ in range(mult):
            addH(pred2, HQ, shift)
        # S_a/(J(Bx)+J(P)) Hilbert function
        boxgens = [{tuple((2 * radii[X] + 1) if a == b else 0 for b in range(D)): 1} for a, X in enumerate(keys)]
        top_deg = sum(2 * radii[X] for X in keys)
        H3 = []
        for i in range(top_deg + 1):
            nmon = len(monos(D, i))
            H3.append(nmon - ideal_dim_in_degree(i, D, boxgens + JP))
        while H3 and H3[-1] == 0:
            H3.pop()
        for _ in range(mult):
            addH(pred3, H3, shift)
        if H3 != HQ:
            bip_fail.append((radii, HQ, H3))
    print('cuts', cuts, 'L', Ls, '|A|', len(A), 'atoms', {X: full[X] for X in keys})
    print('  H_A        ', HA)
    print('  fibred (P2)', pred2, 'equal' if pred2 == HA else 'DIFF')
    print('  R/(x^3+J(P)) (P3)', pred3, 'equal' if pred3 == HA else 'DIFF')
    print('  |P| =', len(P), 'H_P =', HP, ' BIP failures:', len(bip_fail))
    for f in bip_fail[:4]:
        print('    fail', f)
    return HA, pred2, pred3

if __name__ == '__main__':
    tests = [
        (5, [{0, 1, 2}, {2, 3, 4}], [1, 1]),
        (6, [{0, 1, 2, 3}, {2, 3, 4, 5}], [1, 1]),
        (6, [{0, 1, 2, 3}, {2, 3, 4, 5}], [2, 1]),
        (6, [{0, 1, 2}, {1, 2, 3, 4}], [1, 2]),
        (7, [{0, 1, 2, 3}, {3, 4, 5, 6}], [1, 1]),
        (6, [{0, 1, 2, 3}, {2, 3, 4, 5}, {0, 2, 4}], [1, 1, 1]),
        (6, [{0, 1, 2, 3}, {0, 1, 4, 5}, {2, 3, 4, 5}], [1, 1, 1]),
        (7, [{0, 1, 2, 3}, {0, 1, 4, 5}, {0, 2, 4, 6}], [1, 1, 1]),
    ]
    for d, cuts, Ls in tests:
        run(d, cuts, Ls)
        sys.stdout.flush()
