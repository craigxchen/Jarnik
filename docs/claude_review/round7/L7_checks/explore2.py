import time
from fractions import Fraction as Fr
from rfun import RF
import qpoly as Q
import cub13rf as C
x = RF.x()
base = dict(tau=Fr(-2), s34=Fr(3), s35=Fr(5, 2), s45=Fr(7, 3))
for var in ['s45', 's35', 's34', 'tau']:
    t0 = time.time()
    args = dict(base)
    args[var] = x
    e = C.E(args['tau'], args['s34'], args['s35'], args['s45'])
    num = e.n
    rr = Q.rational_roots(num)
    print(var, 'num deg', Q.deg(num), 'den deg', Q.deg(e.d), 'rational roots of num', sorted(rr),
          'time %.1f' % (time.time() - t0), flush=True)
    print('   den rational roots', sorted(Q.rational_roots(e.d)), flush=True)
