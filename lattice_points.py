### EDIT TO COMPUTE ORDER RATHER THAN COUNT

import math
import gmpy2
from gmpy2 import mpz

beta = 1/3
n = 4

# R^2 with 4 lattice points: 548169245, 567454025 

MIN_R = 5
MAX_R = 5

min_coeffs = {}
for rad_squared in range(MIN_R, MAX_R + 1):
    
    radius = math.sqrt(rad_squared)

    # find lattice points
    lp = []
    
    for i in reversed(range(math.ceil(radius/math.sqrt(2)), math.floor(radius)+1)):
        cand = mpz(rad_squared - i**2)
        if cand == 0:
            lp.append((0,radius))
        # for small radius, can compare against list of perfect squares instead?
        elif gmpy2.is_square(cand):
            lp.append((int(math.sqrt(cand)),i))
    
    if not lp:
        continue
    
    for pt in reversed(lp):
        lp.append(pt[::-1])
        
    if radius.is_integer():
        lp.extend([(x,-y) for (x,y) in lp[-2::-1]])
    else:
        lp.extend([(x,-y) for (x,y) in lp[::-1]])
    
    min_coeff = (math.atan2(*lp[0][::-1])-math.atan2(*lp[n-1][::-1])) * rad_squared**(0.5*(1-beta))
    for i in range(1,len(lp)-n-1):
        coeff = (math.atan2(*lp[i][::-1])-math.atan2(*lp[i+n-1][::-1])) * rad_squared**(0.5*(1-beta))
        if min_coeff > coeff:
            min_coeff = coeff
        else:
            continue
    
    min_coeffs[rad_squared] = min_coeff

print(f"minimum ratio for n={n} found: {min(min_coeffs.values())}")