import time, numpy as np
import modla
acc={}
def timed(name,f):
    def g(*a,**k):
        t=time.time(); r=f(*a,**k); acc[name]=acc.get(name,0)+time.time()-t; return r
    return g
modla._split=timed('split',modla._split)
modla._combine=timed('combine',modla._combine)
modla._naive_chunk=timed('naive',modla._naive_chunk)
modla.mm=timed('mm',modla.mm)
orig_reduce=modla.Echelon._reduce
modla.Echelon._reduce=timed('reduce_total',orig_reduce)
modla.Echelon._merge_two=timed('merge_total',modla.Echelon._merge_two)
from combi import *
import iso
orbs=shell_orbits(8,10,6)
T=time.time()
r=iso.lam_profile((4,3,1),orbs,modla.P1)
print('total',time.time()-T, acc)
