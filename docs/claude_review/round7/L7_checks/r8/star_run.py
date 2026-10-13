import sys, random, itertools
from star import *
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
rnd = lambda: R(random.randint(-9, 9), random.randint(1, 5))
planes = star_arrangement(rnd)
# 8 points P_{ai bj ck}: labels a:1,2 b:3,4 c:5,6
trip8 = [(i, j, k) for i in (1,2) for j in (3,4) for k in (5,6)]
chosen = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
pts = [point_of([planes[i], planes[j], planes[k]]) for (i,j,k) in chosen]
x = []
while len(x) < 7:
    c = rnd()
    if c not in x: x.append(c)
data, Q, Lx = member_through(planes, pts, x)
ev, inf, P = analyze(data)
print('contacts at infinity:', inf)
tot = {}
for T, e in sorted(ev.items(), key=lambda kv: (len(kv[0]), kv[0])):
    print([LAB[a] for a in T], 'deg', e.degree(), 'squarefree', e.sqf_part().degree() == e.degree())
pair = sum(e.degree() for T, e in ev.items() if len(T) == 2)
trip = sum(e.degree() for T, e in ev.items() if len(T) in (3, 4))
print('pair-type', pair, 'triple-type', trip, '(B+8K) =', trip - sp.Rational(5,3)*pair, '-K =', sp.Rational(pair,3))
# all event polys pairwise coprime and squarefree:
allp = list(ev.values()); prod = Poly(1, t, domain=QQ)
for e in allp: prod = prod*e
print('total contact points', prod.degree(), 'squarefree product:', prod.sqf_part().degree() == prod.degree())
