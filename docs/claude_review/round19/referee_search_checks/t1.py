import sys; sys.path.insert(0,'.')
from rcore import *
print(full_gain([(0,1)],[(0,2)],1,1,True))
print(full_gain([(0,1)],[(0,2)],-1,1,True))
print(full_gain([(0,1)]*2,[(0,2),(1,2)],-1,4,True))
print(full_gain([(0,1)]*4,[(0,2)]*3+[(1,2)],-27,256,True))
print(full_gain([(0,1)]*4,[(0,2)]*2+[(1,2)]*2,1,16,True))
print(full_gain([(0,1)]*2+[(2,3)]*2,[(0,2),(0,3),(1,2),(1,3)],-1,4,True))
print(full_gain([(0,1)],[(2,3)],1,1,True))
print(full_gain([(0,1)]*5,[(0,4)]+[(1,4)]*4,256,3125,True))
# disc test
p=[3,-7,11,2,-5,13,1]; q=[1,4,-2,9,6,-3,8]
a=restrict([(0,1)]*2,p,q); c=restrict([(0,2),(1,2)],p,q)
print(a,c,disc_poly(a,c,2), rational_roots(disc_poly(a,c,2)))
