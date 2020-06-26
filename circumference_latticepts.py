import math
import gmpy2

epsilon = 0.01

radius = 25
max_ang = 16**(1/3) * radius**(-0.5+epsilon)

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
# angles[0:0] = angles[::-1]

count = []
for i in range(len(lp)):
    # check how many lattice points would fit in arc of length radius^0.5
    running_sum = 0
    running_count = 1

    for ang in angles[i:]:
        running_sum += ang
        if running_sum > max_ang:
            break
        running_count += 1
         
    count.append(running_count)
        