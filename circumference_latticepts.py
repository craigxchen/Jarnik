import math
import gmpy2

radius = 10
max_ang = 1/math.sqrt(radius)

# find lattice points

lp = [(0,radius)]

for i in reversed(range(int(radius//math.sqrt(2)),radius)):
    cand = radius**2 - i**2
    if gmpy2.is_square(cand):
        lp.append((math.sqrt(cand),i))
        # print(i)

angles = []
for j in range(len(lp)-1):
    angles.append( math.atan2(*lp[j][::-1])-math.atan2(*lp[j+1][::-1]) )
    
angles.append( math.atan2(*lp[-1][::-1]) - math.pi/4 )

