import math
import gmpy2

radius = 1000000

# find lattice points

lp = [(0,radius)]

for i in range(int(radius//math.sqrt(2)),radius):
    cand = radius**2 - i**2
    if gmpy2.is_square(cand):
        lp.append((math.sqrt(cand),i))
        # print(i)
    