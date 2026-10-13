# Check arrangement claims of L7 Sec 3.1 and Theorem B class, independently.
import sys; sys.path.insert(0, '.')
from reflib import *
import sympy
# (A1): every six planes independent (rank 4) -> no point on six planes
ok = all(sympy.Matrix([PL[i] for i in S]).rank() == 4 for S in combinations(range(7), 6))
print('(A1) every 6 planes rank 4:', ok)
# coaxial triples: rank 2
coax = [S for S in combinations(range(7), 3) if sympy.Matrix([PL[i] for i in S]).rank() == 2]
print('coaxial triples:', [[LAB[i] for i in S] for S in coax])
# points on five planes: 5-subsets of rank 3
five = [S for S in combinations(range(7), 5) if sympy.Matrix([PL[i] for i in S]).rank() == 3]
print('5-subsets of rank 3 (points on 5 planes):', [[LAB[i] for i in S] for S in five])
four = [S for S in combinations(range(7), 4) if sympy.Matrix([PL[i] for i in S]).rank() == 3]
print('#4-subsets with a common point:', len(four))
# triple points that are not on coaxial triples, i.e. distinct points H_a∩H_b∩H_c
pts = {}
for S in combinations(range(7), 3):
    if sympy.Matrix([PL[i] for i in S]).rank() == 3:
        x = tuple(meet3(*S))
        # normalise sign
        for c in x:
            if c != 0:
                if c < 0: x = tuple(-y for y in x)
                break
        pts.setdefault(x, []).append(S)
print('#distinct triple points:', len(pts))
for x, Ss in pts.items():
    planes = set(i for S in Ss for i in S)
    if len(planes) == 3:
        print('  simple triple point', x, [LAB[i] for i in sorted(planes)])
cube = {}
for i in (1, 2):
    for j in (3, 4):
        for k in (5, 6):
            cube[(i, j, k)] = meet3(i, j, k)
Ekeys = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
E = [cube[k] for k in Ekeys]
webm = web(E)
print('dim of quadrics through E:', len(webm))
