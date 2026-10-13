import sys, random, itertools
from star import *
random.seed(2)
rnd = lambda: R(random.randint(-9, 9), random.randint(1, 5))
planes = star_arrangement(rnd)
trip8 = [(i, j, k) for i in (1,2) for j in (3,4) for k in (5,6)]
pts8 = {(i,j,k): point_of([planes[i], planes[j], planes[k]]) for (i,j,k) in trip8}
x = [rnd() for a in range(7)]
v = line_through_x(planes, x)[0]; L1 = v[:4,:]; L2 = v[4:,:]
for S in itertools.combinations(trip8, 6):
    Qs = quadric_through([L1, L2, L1+L2] + [pts8[s] for s in S])
    print(S, 'num', len(Qs), 'ranks', [Q.rank() for Q in Qs])
