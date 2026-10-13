# Probe Sing(F_A) at interior points: the lines m = P1 ∩ P2 for partitions of E into two triples
# (rank-2 web quadric P1 ∪ P2). Test via Lemma 5.3: comp on good mu-curves through the point vs generic -10.
import sys, random; sys.path.insert(0, '.')
from reflib import *
import sympy as sp
cube = {(i,j,k): meet3(i,j,k) for i in (1,2) for j in (3,4) for k in (5,6)}
Ekeys = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
Ep = [cube[k] for k in Ekeys]
webm = web(Ep)
def nm(k): return ''.join(LAB[i] for i in k)
rng = random.Random(5)
for T1 in combinations(range(6), 3):
    if 0 not in T1: continue
    T2 = tuple(i for i in range(6) if i not in T1)
    P1 = rref_null([Ep[i] for i in T1], 4); P2 = rref_null([Ep[i] for i in T2], 4)
    if len(P1) != 1 or len(P2) != 1:
        print('collinear triple', T1); continue
    h1 = intvec(P1[0]); h2 = intvec(P2[0])
    arr1 = [LAB[a] for a in range(7) if sp.Matrix([h1, list(PL[a])]).rank() == 1]
    arr2 = [LAB[a] for a in range(7) if sp.Matrix([h2, list(PL[a])]).rank() == 1]
    mline = rref_null([h1, h2], 4)
    p = intvec(mline[0]); q = intvec(mline[1])
    pos = [(ev(PL[a], q), -ev(PL[a], p)) for a in range(7)]
    coll = [(LAB[i], LAB[j]) for i, j in combinations(range(7), 2) if pos[i][0] * pos[j][1] - pos[i][1] * pos[j][0] == 0]
    inpl = [LAB[a] for a in range(7) if pos[a] == (0, 0)]
    # is the web quadric through m unique, and what is it
    rows = [[sum(Fr(S[i][j]) * p[i] * p[j] for i in range(4) for j in range(4)) for S in webm],
            [sum(Fr(S[i][j]) * q[i] * q[j] for i in range(4) for j in range(4)) for S in webm],
            [sum(Fr(S[i][j]) * p[i] * q[j] for i in range(4) for j in range(4)) for S in webm]]
    nsp = rref_null(rows, 4)
    line = 'T1=%s T2=%s planes %s|%s; m: collisions %s inplanes %s; #web quadrics through m: %d' % (
        [nm(Ekeys[i]) for i in T1], [nm(Ekeys[i]) for i in T2], arr1, arr2, coll, inpl, len(nsp))
    print(line)
    if coll or inpl:
        continue
    # interior point x_m: positions as affine numbers
    xm = [Fr(a, b) if b != 0 else None for (a, b) in pos]
    # good mu-curves through x_m: V_a(s) = (x_a*den + s*r1, den + s*r2) style: g_a(0) = x_m,a
    vals = []
    tries = 0
    while len(vals) < 3 and tries < 200:
        tries += 1
        V = []
        for (A, Bb) in pos:
            A = Fr(A); Bb = Fr(Bb)
            L = A.denominator * Bb.denominator
            V.append([P([int(A * L), rng.randint(-9, 9)]), P([int(Bb * L), rng.randint(-9, 9)])])
        cls, good, info = class_and_goodness(V)
        if not good or not all(cls.get(frozenset(k), 0) == 2 for k in combinations(range(7), 2)) or sum(cls.values()) != 42:
            continue
        Br, dBr, _ = FA_brackets(V, webm)
        c, nk = comp_count(Br, dBr)
        vals.append(c)
    print('    comp_A on good mu-curves through x_m:', vals, '(generic value -10)')
