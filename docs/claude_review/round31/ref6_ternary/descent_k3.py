"""Direct product over F_729[t] for q=27, k=3 (19683 linear factors) at small L, compared with marks.py;
checks the F_9 descent and the norm identity F*conj(F) = P_V^2 for k=3."""
import numpy as np, itertools, sys, time
import descent_check as dc
from marks import raw_mark
from f9lib import F9P
K = dc.K
L = int(sys.argv[1])
forms = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1)]
signs = [1,1,1,1,1]
PV = None
for lam, e in zip(forms, signs):
    t0 = time.time()
    p = dc.product_point(lam, e, L, k=3)
    z = dc.to_f9(p)
    assert z is not None
    F, A, c = raw_mark(lam, e, 27, F9P([0]*L + [1]))
    assert F == z
    nz = z * z.conj()
    if PV is None: PV = nz
    assert nz == PV and nz.is_f3()
    print(lam, 'direct product == closed form; coefficients in F_9; norm equal, in F_3[t]; deg', z.deg(), '%.1fs' % (time.time()-t0), flush=True)
