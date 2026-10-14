"""Lattice structure of Sylvester column-family images.

t-dim XOR block U = F_2^t (characters a != 0), classes z in F_2^t with class sums c_z (sum c = 0).
Column a: y_a = k(S_a) = sum_{z: a.z = 1} c_z  (the sign of the convention y = -k(S_a) is irrelevant).
Lambda_t = image of Z^{2^t}_0;  for F subset U\\0: Lambda_F = projection to coordinates F.
Reports: Smith invariants of Lambda_F in Z^F (index), whether Lambda_F is a 'binary code lattice'
(2 Z^F subset Lambda_F), the code C_F = Lambda_F^perp mod 2 and its minimum distance, and (small j)
whether proj_F(Y_t(j)) = Lambda_F cap [-j,j]^F.
"""
import itertools, sys, random
from flint import fmpz_mat

def dot(a, z):
    return bin(a & z).count('1') % 2

def gens(t, F):
    z0 = 0
    rows = []
    for z in range(1, 2 ** t):
        rows.append([dot(a, z) - dot(a, z0) for a in F])
    return rows

def snf_diag(rows, b):
    M = fmpz_mat(rows)
    # Smith form via flint: use hnf then compute invariants by gcd of minors is heavy; use sympy for small
    from sympy import Matrix
    from sympy.matrices.normalforms import smith_normal_form
    S = smith_normal_form(Matrix(rows))
    diag = [abs(S[i, i]) for i in range(min(S.shape)) if S[i, i] != 0]
    return diag

def code_of(t, F):
    """binary code C = {s in F_2^F : <s, y> even for all y in Lambda_F}"""
    rows = gens(t, F)
    b = len(F)
    C = []
    for s in range(1, 2 ** b):
        if all(sum(r[q] for q in range(b) if (s >> q) & 1) % 2 == 0 for r in rows):
            C.append(s)
    return C

def lattice_member(t, F, y):
    """y in Lambda_F ?  brute force over extensions is expensive; use rank test via integer solve"""
    from sympy import Matrix
    raise NotImplementedError

def image_equals_section(t, F, j):
    """compare proj_F(Y_t(j)) with Lambda_F cap box, by enumeration (small only)."""
    allc = list(range(1, 2 ** t))
    mod = 2 ** (t - 1)
    Y = set()
    for y in itertools.product(range(-j, j + 1), repeat=len(allc)):
        ok = True
        for z in range(2 ** t):
            s = sum(yy * (1 - 2 * dot(a, z)) for yy, a in zip(y, allc))
            if s % mod:
                ok = False; break
        if ok:
            Y.add(tuple(y[allc.index(a)] for a in F))
    # Lambda_F cap box: Lambda_F = proj of Lambda_t; Lambda_t = {y : sum y_a chi_a(z) = 0 mod 2^{t-1}}
    # membership in proj: exists extension in Z (not boxed): test with extension range wide enough
    rest = [a for a in allc if a not in F]
    R = 2 ** (t - 1)
    sec = set()
    for yF in itertools.product(range(-j, j + 1), repeat=len(F)):
        found = False
        for yr in itertools.product(range(0, R), repeat=len(rest)):  # residues suffice (2^{t-1} Z^b in Lambda_t)
            full = dict(zip(F, yF)); full.update(zip(rest, yr))
            if all(sum(full[a] * (1 - 2 * dot(a, z)) for a in allc) % mod == 0 for z in range(2 ** t)):
                found = True; break
        if found:
            sec.add(yF)
    return Y == sec, len(Y), len(sec)

if __name__ == '__main__':
    for t in (2, 3, 4):
        allc = list(range(1, 2 ** t))
        seen = {}
        # orbit representatives up to GL(t,2) would be ideal; sample subsets by size, dedupe by invariants
        for b in range(2, len(allc) + 1):
            subs = list(itertools.combinations(allc, b))
            if len(subs) > 400:
                random.seed(b); subs = random.sample(subs, 400)
            stats = {}
            for F in subs:
                rows = gens(t, list(F))
                diag = snf_diag(rows, b)
                full_rank = len(diag) == b
                C = code_of(t, list(F))
                binary = full_rank and all(dd in (1, 2) for dd in diag)
                dmin = min((bin(s).count('1') for s in C), default=None)
                key = (tuple(diag), binary, len(C), dmin)
                stats[key] = stats.get(key, 0) + 1
            for key, cnt in sorted(stats.items()):
                print(f't={t} b={b}: SNF={list(key[0])} binary={key[1]} |C|={key[2]} dmin={key[3]}  (#{cnt})')
        sys.stdout.flush()
    # image vs section for t=3, small j
    for j in (1, 2):
        for F in [(1, 2, 3), (1, 2, 4), (1, 2, 3, 4), (1, 2, 4, 7), (1, 2, 3, 4, 5), (1, 2, 3, 4, 5, 6), (1, 2, 3, 4, 5, 6, 7)]:
            print('t=3 j=%d F=%s image==Lambda_F cap box:' % (j, F), image_equals_section(3, list(F), j))
            sys.stdout.flush()
