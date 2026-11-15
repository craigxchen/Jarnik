import sys; sys.path.insert(0,'.')
from localgcd import *
from cube_basic import cube
y=[Fr(-9,4),Fr(-2,3),Fr(1,4),Fr(8,3)]
roots,P=cube(y)
pr=[[roots[k].conj() if v[k] else roots[k] for k in range(8)] for (z,v) in P.values()]
print(C_inf(pr,8,verbose=True))
