# Control for t11: same limit test at a generic line m' meeting l_2 (generic point of D_{O b1 b2}) and at a generic line.
import sys, random; sys.path.insert(0, '.')
exec(open('t11_dirlimit.py').read().split('rng = random.Random(3)')[0])
rng = random.Random(4)
l2 = rref_null([PL[0], PL[3]], 4)     # axis l_2 = H_O ∩ H_b1
pt = [Fr(3) * l2[0][i] + Fr(-2) * l2[1][i] for i in range(4)]; pt = intvec(pt)
for name, (A, Bv) in {'generic line meeting l_2': (pt, [5, -3, 7, 2]), 'generic line': ([1, 4, -2, 3], [5, -3, 7, 2])}.items():
    m0, m1 = A, Bv
    res = []
    for t in range(5):
        al, be = rng.randint(1, 9), rng.randint(-9, 9); ga, de = rng.randint(-9, 9), rng.randint(1, 9)
        if al * de - be * ga == 0: continue
        y = [al * m0[i] + be * m1[i] for i in range(4)]; z0 = [ga * m0[i] + de * m1[i] for i in range(4)]
        w = [rng.randint(-20, 20) for _ in range(4)]
        res.append(limdir(y, z0, w))
    print(name, ': distinct limit directions', len(set(res)), 'out of', len(res))
