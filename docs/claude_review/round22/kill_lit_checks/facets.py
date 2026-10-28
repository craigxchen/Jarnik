# Root polytope P = conv(e_i - e_j) in M = {sum=0} (k=5): every 1_A (A proper nonempty) is a facet normal
# (face <u,1_A> = 1 has affine dimension 3 = dim P - 1), and P is reflexive (all facets at distance 1).
# Hence the toric variety X' of the normal fan of P has the same 30 rays as the permutohedral fan; X -> X' is a small
# crepant toric morphism, -K_X = pullback of the ample Cartier -K_{X'}; X' is a (twisted) Gorenstein toric Fano 4-fold.
import itertools, sympy as sp
k = 5
V = [tuple((1 if t == i else 0) - (1 if t == j else 0) for t in range(k)) for i in range(k) for j in range(k) if i != j]
cnt = 0
for r in range(1, k):
    for A in itertools.combinations(range(k), r):
        face = [v for v in V if sum(v[i] for i in A) == 1]
        assert all(sum(v[i] for i in A) <= 1 for v in V)
        M = sp.Matrix([[a - b for a, b in zip(v, face[0])] for v in face[1:]])
        if M.rank() == 3: cnt += 1
print("facets among the 30 inequalities <u,1_A> <= 1:", cnt, "of 30")
# simple polytope? count facets through each vertex (simple <=> exactly 4)
deg = [sum(1 for r in range(1, k) for A in itertools.combinations(range(k), r) if sum(v[i] for i in A) == 1) for v in V]
print("facets through each vertex:", sorted(set(deg)), "(simple iff all = 4)")
