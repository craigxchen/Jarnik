from faithful_search import *
import time, itertools
def cycle(d): return [[i, (i+1) % d] for i in range(d)]
petersen = [[i, (i+1)%5] for i in range(5)] + [[5+i, 5+(i+2)%5] for i in range(5)] + [[i, i+5] for i in range(5)]
cube3 = [[a, b] for a in range(8) for b in range(8) if a < b and bin(a ^ b).count('1') == 1]
K33 = [[a, b] for a in range(3) for b in range(3, 6)]
K44 = [[a, b] for a in range(4) for b in range(4, 8)]
cases = [('C%d' % d, d, cycle(d)) for d in range(5, 10)] + [('K33', 6, K33), ('cube3', 8, cube3), ('K44', 8, K44)]
if len(sys.argv) > 1: cases.append(('petersen', 10, petersen))
for name, d, cols in cases:
    t = time.time()
    r = analyze(d, cols, [1]*len(cols), extra_funcs=True)
    print(name, r, '%.1fs' % (time.time()-t), flush=True)
