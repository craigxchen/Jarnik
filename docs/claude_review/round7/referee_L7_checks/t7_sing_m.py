# Is x_m (m = P1∩P2, P1=<a1b1c1,a1b2c2,a2b2c1>, P2=<a1b1c2,a1b2c1,a2b1c1>) a singular point of F_A?
# m meets l_2 (O,b1,b2 collide) and nothing else: x_m lies on D_{O b1 b2} only (not in ∂^(2)).
# Test (Lemma 5.3): pencils of lines (y,Π) are Mobius curves. Compare comp_A for pencils through m with generic pencils
# of the SAME class (class computed exactly).
import sys, random; sys.path.insert(0, '.')
from reflib import *
cube = {(i,j,k): meet3(i,j,k) for i in (1,2) for j in (3,4) for k in (5,6)}
Ekeys = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
Ep = [cube[k] for k in Ekeys]
webm = web(Ep)
T1, T2 = (0, 3, 5), (1, 2, 4)
h1 = intvec(rref_null([Ep[i] for i in T1], 4)[0]); h2 = intvec(rref_null([Ep[i] for i in T2], 4)[0])
mb = rref_null([h1, h2], 4); m0 = intvec(mb[0]); m1 = intvec(mb[1])
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
def pencil(y, z0, z1):
    zs = [P([z0[i], z1[i]]) for i in range(4)]
    yy = [P([y[i]]) for i in range(4)]
    return [[lin(PL[a], zs), -lin(PL[a], yy)] for a in range(7)]
def rnd(): return [rng.randint(-20, 20) for _ in range(4)]
def keyc(cls): return tuple(sorted(((tuple(sorted(S)), d) for S, d in cls.items())))
res = {}
for trial in range(12):
    # pencil through m: y on m, plane Π = span(m, w); z0 = other point of m, z1 = w  (so s=0 gives R=m)
    al, be = rng.randint(1, 9), rng.randint(-9, 9)
    y = [al * m0[i] + be * m1[i] for i in range(4)]
    ga, de = rng.randint(-9, 9), rng.randint(1, 9)
    z0 = [ga * m0[i] + de * m1[i] for i in range(4)]
    w = rnd()
    V = pencil(y, z0, w)
    cls, good, info = class_and_goodness(V)
    if not good: continue
    try:
        Br, dBr, _ = FA_brackets(V, webm)
    except AssertionError as ex:
        print('degenerate tangent system on pencil through m: y=', y, 'z0=', z0, 'w=', w); continue
    c, nk = comp_count(Br, dBr)
    res.setdefault(keyc(cls), {'thr_m': [], 'generic': []})['thr_m'].append(c)
    # generic pencil: random y, random plane; same contact structure checked by class key
    for tr2 in range(3):
        V2 = pencil(rnd(), rnd(), rnd())
        cls2, good2, _ = class_and_goodness(V2)
        if not good2: continue
        try:
            Br2, dBr2, _ = FA_brackets(V2, webm)
        except AssertionError:
            print('degenerate on generic pencil'); continue
        c2, _ = comp_count(Br2, dBr2)
        res.setdefault(keyc(cls2), {'thr_m': [], 'generic': []})['generic'].append(c2)
for k, v in res.items():
    print('class', [('-'.join(LAB[i] for i in S), d) for S, d in k])
    print('   comp through x_m:', v['thr_m'], '  comp generic pencils:', v['generic'])
