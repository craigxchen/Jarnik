import sys, random, itertools
from star import *
random.seed(2)
rnd = lambda: R(random.randint(-9, 9), random.randint(1, 5))
planes = star_arrangement(rnd)
chosen = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
pts = [point_of([planes[i], planes[j], planes[k]]) for (i,j,k) in chosen]
x = [rnd() for a in range(7)]
data, Q, Lx = member_through(planes, pts, x)
print('Q rank', Q.rank())
for a,(N,D) in enumerate(data): print(LAB[a], N, '/', D)
n=7; P={}
for a, b in itertools.combinations(range(n), 2):
    Na, Da = data[a]; Nb, Db = data[b]
    P[(a,b)] = Poly(sp.expand(Na*Db - Nb*Da), t, domain=QQ)
    print(LAB[a],LAB[b], sp.factor(P[(a,b)].as_expr()))
