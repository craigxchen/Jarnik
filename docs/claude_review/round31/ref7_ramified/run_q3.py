import itertools
from fp import *
from construct_q3 import *
Delta = trim([1,0,2])          # 2t^2+1 = 2(t-1)(t+1)
C = trim([0,2,0,1])            # t^3 - t, divisible by (t-1)(t+1) and by t
# best signs for each m
for m in (2,3,4):
    for forms in itertools.combinations(ALLFORMS, m):
        best=None
        for signs in itertools.product((1,-1), repeat=m):
            o = build(Delta, C, 1, list(forms), list(signs), check_sf=False)
            if best is None or o['g']>best[0]: best=(o['g'],signs)
        print(m, forms, 'best g', best)
