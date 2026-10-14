import time
from flib import *
for name, fn, arg in [('sdc7', sdc_patterns, 7), ('pal7', pal_patterns, 7), ('syl8', sylvester_patterns, 3),
                      ('sdc11', sdc_patterns, 11), ('pal11', pal_patterns, 11)]:
    M, pats = fn(arg)
    for L in (1, 2, 3):
        t = time.time(); A = box(M, pats, L)
        print(name, 'M=%d #pats=%d L=%d |A|=%d (%.1fs)' % (M, len(pats), L, len(A), time.time() - t), flush=True)
