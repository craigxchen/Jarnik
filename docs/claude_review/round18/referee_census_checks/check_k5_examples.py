from rlib import *
from math import log
for x in [(0,22,35,46,85),(0,39,50,63,85),(0,17,29,48,99),(0,51,70,82,99)]:
    d=contacts_classes(list(x))
    lg=[log(v) if v>1 else 0 for v in d.values()]
    lo,hi=min(lg),max(lg)
    eps=(hi-lo)/(hi+lo) if lo>0 else 1
    print(x,"cov",sum(v>1 for v in d.values()),"eps %.4f w* %.3f"%(eps,(hi+lo)/2), sorted(d.values()))
