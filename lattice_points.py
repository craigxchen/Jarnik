import math

beta = 1/2

MIN_R = 5695325
MAX_R = 5695325

max_counts = {}
for radius in range(MIN_R, MAX_R + 1):
    radius = math.sqrt(radius)
    
    arc_length = radius**(beta) 
    max_ang = arc_length / radius
    
    # find lattice points
    lp = []
    
    for i in reversed(range(math.ceil(radius//math.sqrt(2)), math.floor(radius))):
        cand = radius**2 - i**2
        if cand == 0:
            lp.append(0,radius)
        # for small radius, can compare against list of perfect squares instead?
        elif math.sqrt(cand).is_integer():
            lp.append((math.sqrt(cand),i))
    
    # angle from first point to its reflection across y-axis
    angles = [2 * math.atan2(*lp[0])]
    
    # angles between lattice points
    for j in range(len(lp)-1):
        angles.append( math.atan2(*lp[j][::-1])-math.atan2(*lp[j+1][::-1]) )
        
    # angle between point closest to y=x and its reflection
    angles.append( 2*math.atan2(*lp[-1][::-1]) - math.pi/2 )
    # extend to first quadrant
    angles.extend(angles[-2::-1])
    # extend to entire right semi-circle
    angles.extend(angles)

    
    count = []
    for i in range(len(angles)):
        # check how many lattice points would fit in arc of length radius*max_ang
        running_sum = 0
        running_count = 1
    
        for ang in angles[i:]:
            running_sum += ang
            if running_sum >= max_ang:
                break
            running_count += 1
             
        count.append(running_count)
    
    max_counts[radius] = max(count)

print(f"maximum lattice points found: {max(max_counts.values())}")