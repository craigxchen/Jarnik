from sympy import Matrix, eye, zeros
import itertools
edges=[(i,j) for i in range(5) for j in range(i+1,5)]
M=Matrix(10,10,lambda a,b: 1 if a!=b and len(set(edges[a])&set(edges[b]))==1 else 0)
print("det line-graph adjacency K5:",M.det(), " 6*M^-1 integral:", all(x.is_integer for x in 6*M.inv()))
v=[Matrix([1 if k==i else 0 for k in range(4)]) for i in range(4)]+[Matrix([-1]*4)]
b=[v[i]+v[j] for (i,j) in edges]
mons=[(a,c) for a in range(4) for c in range(a,4)]
Q=Matrix(10,10,lambda r,s: b[r][mons[s][0]]*b[r][mons[s][1]])
print("quadric evaluation det:",Q.det())
# B diag(x) Z rank-drop structure: B = 4x10 pair-column matrix; kernel dim
B=Matrix.hstack(*b); print("rank B",B.rank(), "kernel dim", len(B.nullspace()))
