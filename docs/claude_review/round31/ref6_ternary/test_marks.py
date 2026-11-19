import pickle, itertools
from marks import *
from f9lib import gcd9, divmod9
d = pickle.load(open('four_mark_L11.pkl','rb'))
forms = [(1,0),(0,1),(1,1),(1,2)]
M = check_family(forms, [1]*4, 27, 11)
G = M[0]
for z in M[1:]: G = gcd9(G, z)
red = [divmod9(z, G)[0] for z in M]
from f9lib import F9P
assert all(F9P(*a) == b for a, b in zip(d['red'], red))
print('general-k subspace-polynomial code reproduces the direct-product four-mark points')
