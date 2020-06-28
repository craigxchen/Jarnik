import math
import gmpy2

epsilon = 0

MIN_R = 1000
MAX_R = 10000

max_counts = {}
for radius in range(MIN_R, MAX_R + 1):
    
    arc_length = math.sqrt(2) * radius**(0.5-epsilon)
    max_ang = arc_length / radius
    
    # find lattice points
    lp = [(0,radius)]
    
    for i in reversed(range(int(radius//math.sqrt(2))+1,radius)):
        cand = radius**2 - i**2
        # for small radius, can compare against list of perfect squares instead?
        if gmpy2.is_square(cand):
            lp.append((math.sqrt(cand),i))
    
    angles = []
    for j in range(len(lp)-1):
        angles.append( math.atan2(*lp[j][::-1])-math.atan2(*lp[j+1][::-1]) )
        
    angles.append( 2*math.atan2(*lp[-1][::-1]) - math.pi/2 )
    angles.extend(angles[-2::-1])
    angles.extend(angles)

    
    count = []
    for i in range(len(angles)):
        # check how many lattice points would fit in arc of length radius*max_ang
        running_sum = 0
        running_count = 1
    
        for ang in angles[i:]:
            running_sum += ang
            if running_sum > max_ang:
                break
            running_count += 1
             
        count.append(running_count)
    
    max_counts[radius] = max(count)

print(f"maximum lattice points found: {max(max_counts.values())}")