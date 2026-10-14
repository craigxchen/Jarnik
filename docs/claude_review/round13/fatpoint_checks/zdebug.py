import itertools, numpy as np
from fplib import *
M, pats = MODELS['sdc7'](); A = box(M, pats, 2); d = 7; m = 6
G = [(0, 1, 0, 0, 0, 0, -1), (0, 0, 0, 0, 1, -1, 0), (0, 1, 0, 0, -1, 0, 0), (1, 0, 0, 0, 0, -1, 0), (1, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 1)]
S = np.array(G); Aset = {a:i for i,a in enumerate(A)}
verts = []
for r in range(m+1):
    for I in itertools.combinations(range(m), r):
        verts.append((S[list(I)].sum(axis=0) if I else np.zeros(d,dtype=int), (-1)**r))
lo = np.min([v for v,_ in verts],axis=0); hi = np.max([v for v,_ in verts],axis=0)
for k0 in itertools.product(*[range(-1-int(l), 2-int(h)) for l,h in zip(lo,hi)]):
    if all(tuple(v+np.array(k0)) in Aset for v,_ in verts): break
c = np.zeros(len(A), dtype=np.int64)
for v,s in verts: c[Aset[tuple(v+np.array(k0))]] += s
low = []
for e in range(m): low += reduced_monomials(d,e)
print('E_low c == 0:', not (eval_matrix(A, low) @ c).any())
monos_m = reduced_monomials(d, m)
mu = eval_matrix(A, monos_m) @ c
poly = product_coeffs([list(g) for g in G], m)
ratio = set()
for a,i in zip(monos_m, range(len(monos_m))):
    pc = poly.get(a, 0)
    if mu[i] or pc:
        ratio.add((mu[i]*multinom(m,a)) / pc if pc else 'inf')
print('ratios mu*multinom/p:', ratio)
nonred = [e for e in poly if max(e)>2]
print('nonreduced exps in product:', nonred[:5])
monos_m2, Tb = tangent_space_mod(A, m)
print('dimT', Tb.shape)
rows = Tb.tolist() + [[int(x) % PR for x in mu]]
print('mu in T mod p:', nmod_mat(len(rows), len(monos_m), [x for r in rows for x in r], PR).rank() == Tb.shape[0])
