import sys, random, itertools
from star import *
random.seed(1)
rnd = lambda: R(random.randint(-9, 9), random.randint(1, 5))
planes = star_arrangement(rnd)
trip8 = [(i, j, k) for i in (1,2) for j in (3,4) for k in (5,6)]
pts8 = [point_of([planes[i], planes[j], planes[k]]) for (i,j,k) in trip8]
M8 = Matrix.hstack(*pts8)
print('rank of 8 points', M8.rank())
for S in itertools.combinations(range(8), 4):
    if Matrix.hstack(*[pts8[s] for s in S]).rank() < 4: print('coplanar 4:', [trip8[s] for s in S])
# quadrics through all 8 points
print('quadrics through 8 pts:', len(quadric_through(pts8)))
print('quadrics through first 6:', len(quadric_through(pts8[:6])))
x = [rnd() for a in range(7)]
v = line_through_x(planes, x)[0]; L1 = v[:4,:]; L2 = v[4:,:]
Qs = quadric_through([L1, L2, L1+L2] + pts8[:6]); Q = Qs[0]
print('rank Q', Q.rank(), Q)
