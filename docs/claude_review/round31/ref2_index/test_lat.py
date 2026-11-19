import random
from sympy import Matrix
from lat import lattice_index
random.seed(1)
bad=0
for t in range(400):
    k=random.randint(1,4); n=random.randint(k,7)
    vecs=[[random.randint(-6,6) for _ in range(k)] for _ in range(n)]
    A=Matrix(vecs).T  # k x n
    r=A.rank()
    if r<k: ref=0
    else:
        from itertools import combinations
        from math import gcd
        g=0
        for cs in combinations(range(n),k):
            g=gcd(g,int(abs(A[:,list(cs)].det())))
        ref=g
    if lattice_index(vecs,k)!=ref: bad+=1
print("mismatches",bad)
print(lattice_index([[7,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],4))
