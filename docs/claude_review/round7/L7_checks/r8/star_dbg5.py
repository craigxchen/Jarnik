import sys, random, itertools
from star import *
random.seed(2)
rnd = lambda: R(random.randint(-9, 9), random.randint(1, 5))
planes = star_arrangement(rnd)
chosen = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
pts = [point_of([planes[i], planes[j], planes[k]]) for (i,j,k) in chosen]
x = [rnd() for a in range(7)]
print('x', x)
v = line_through_x(planes, x)[0]; L1 = v[:4,:]; L2 = v[4:,:]
for a in range(7):
    pa = x[a]*L1 + L2
    print(LAB[a], 'on plane', (planes[a]*pa)[0])
Q = quadric_through([L1, L2, L1+L2] + pts)[0]
lab = kernel(Matrix.vstack(planes[1], planes[4]))
print('axis a1b2 points on Q:', [(k.T*Q*k)[0] for k in lab], (lab[0].T*Q*lab[1])[0])
# does R_x meet axis a1 b2?
print('R_x meets axis a1b2:', Matrix.hstack(L1, L2, lab[0], lab[1]).rank())
print('x_a1, x_b2', x[1], x[4])
