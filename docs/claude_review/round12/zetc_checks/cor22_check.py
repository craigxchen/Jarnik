# Corollary 2.2(3): hexagon cut S={0,1}, L=1, d=6: rho(A)=11; nu(u)=7 for u with distinct nonzero u_0,u_1, two more nonzero coords
from tclib import box_points, nu_exact
pts = box_points(6, [([0,1], 1)])
for u in [(1,2,3,5,0,0), (2,7,-3,11,0,0), (1,2,1,1,0,0)]:
    print(u, 'nu, rho =', nu_exact(pts, u))
