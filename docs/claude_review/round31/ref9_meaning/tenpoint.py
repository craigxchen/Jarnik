"""Own check of the ten-point illustration of IX.2 and of the orientation reading."""
import math
from sympy import factorint
from itertools import combinations
Z=[(24791,23712),(24745,23760),(24636,23873),(24567,23944),(24260,24255),(24255,24260),
   (23944,24567),(23873,24636),(23760,24745),(23712,24791)]
N=Z[0][0]**2+Z[0][1]**2
assert all(x*x+y*y==N for x,y in Z)
print("N",N,factorint(N))
g=0
for x,y in Z: g=math.gcd(g,math.gcd(x,y))
print("rational gcd of all coords",g)
R=math.sqrt(N); ch=max(math.hypot(a[0]-b[0],a[1]-b[1]) for a,b in combinations(Z,2))
print("chord constant",ch/N**0.25)
# reflection pairs z_j = i*conj(z_i) (swap x,y): T(c)=1 certificates
pairs=[(i,j) for i,j in combinations(range(10),2) if Z[j]==(Z[i][1],Z[i][0])]
print("swap pairs",pairs)
xs=sorted(set(x+y for x,y in Z)); print("projections Re(conj(1+i)z)=x+y:",xs,"diameter",xs[-1]-xs[0])
th=sum(math.atan2(y,x) for x,y in Z)/10
print("midpoint angle - pi/4 =",th-math.pi/4,"; times sqrt(R) =",(th-math.pi/4)*math.sqrt(R))
