"""EXACT: boundary-intersection vectors of curve classes given in the Kapranov psi_7 model
(e; m_a at p_a, m_ab along L_ab, m_abc along Pi_abc), with D_{a7}=E_a, D_{ab7}=E_ab, D_{abc7}=E_abc,
D_{mn} (m,n<=6) = hyperplane through p_I, I=[6]-{m,n}.  Prints pair, triple, -K, B+8K and the ratio
(B+8K)/B for the fibre classes of the known Q-rational fibrations, cub13, cub12, beta.
Also the toric 1-PS classes gamma_{a,b,A} computed directly from contacts."""
import itertools
from fractions import Fraction as Fr
L6 = range(1, 7)
def vec(e, ma={}, mab={}, mabc={}):
    D = {}
    for a in L6: D[frozenset({a, 7})] = ma.get(a, 0)
    for ab in itertools.combinations(L6, 2): D[frozenset(ab) | {7}] = mab.get(ab, 0)
    for abc in itertools.combinations(L6, 3): D[frozenset(abc) | {7}] = mabc.get(abc, 0)
    for mn in itertools.combinations(L6, 2):
        I = [x for x in L6 if x not in mn]
        D[frozenset(mn)] = (e - sum(ma.get(a, 0) for a in I) - sum(mab.get(ab, 0) for ab in itertools.combinations(I, 2))
                            - sum(mabc.get(abc, 0) for abc in itertools.combinations(I, 3)))
    return D
def norm(S):
    S = frozenset(S); C = frozenset(range(1, 8)) - S
    return S if len(S) < len(C) or (len(S) == len(C) and min(S) < min(C)) else C
def stats(D):
    pair = sum(v for k, v in D.items() if len(norm(k)) == 2)
    trip = sum(v for k, v in D.items() if len(norm(k)) == 3)
    return pair, trip, Fr(pair, 3), trip - Fr(5, 3) * pair
def onepS(a, b, A):
    """1-PS orbit in torus T_ab with direction e_A (A subset of the other five labels)"""
    D = {}
    rest = [x for x in range(1, 8) if x not in (a, b)]
    D[norm({a} | set(A))] = D.get(norm({a} | set(A)), 0) + 1
    D[norm({b} | set(A))] = D.get(norm({b} | set(A)), 0) + 1
    for i in A:
        for j in rest:
            if j not in A: D[norm({i, j})] = D.get(norm({i, j}), 0) + 1
    return D
planes13 = [(1,2,3),(1,2,4),(1,2,5),(1,2,6),(1,3,4),(1,3,5),(1,3,6),(1,4,5),(2,4,6),(2,5,6),(3,4,6),(3,5,6),(4,5,6)]
classes = {
  'beta (balanced)': None,
  'f_1 = line through p_1 (forgetful fibre)': vec(1, ma={1: 1}),
  'f_7 = RNC through p_1..p_6 (forgetful fibre)': vec(4, ma={a: 1 for a in L6}),
  'general line (projection from a general point)': vec(1),
  'line meeting 3 planes 123,145,246 (plane projections)': vec(1, mabc={(1,2,3): 1, (1,4,5): 1, (2,4,6): 1}),
  'line meeting L_12 and Pi_345': vec(1, mab={(1,2): 1}, mabc={(3,4,5): 1}),
  'cub12 (Q-rational 4-dim family)': vec(3, mabc={p: 1 for p in planes13 if p != (4,5,6)}),
  'cub13': vec(3, mabc={p: 1 for p in planes13}),
  '1-PS gamma_{6,7,{1}} (toric)': onepS(6, 7, [1]),
  '1-PS gamma_{6,7,{1,2}} (toric)': onepS(6, 7, [1, 2]),
}
for name, D in classes.items():
    if D is None:
        print('%-58s pair 21 triple 35  -K 7  B+8K 0' % name); continue
    p, t, k, b8 = stats(D)
    assert all(v >= 0 for v in D.values()), name
    print('%-58s pair %2d triple %2d  -K %s  B+8K %s  ratio (B+8K)/B = %s' % (name, p, t, k, b8, b8 / (p + t)))
