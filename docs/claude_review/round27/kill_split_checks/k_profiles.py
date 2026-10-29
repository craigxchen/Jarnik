# Pi_nt dimension, one-sided vertices, code hull; own code
import itertools, sympy as sp
from fractions import Fraction as Fr
Bs=[frozenset(s) for r in range(1,5) for s in itertools.combinations(range(5),r)]
idx={B:i for i,B in enumerate(Bs)}; full=frozenset(range(5))
cuts=[B for B in Bs if 0 in B]  # representative containing 0 -> 15 cuts
assert len(Bs)==30 and len(cuts)==15
rows=[]
for B in cuts:
    r=[0]*30; r[idx[B]]=1; r[idx[full-B]]=1; rows.append(r)
for a in range(1,5):
    r=[0]*30
    for B in Bs:
        r[idx[B]]+= (1 if a in B else 0) - (1 if 0 in B else 0)
    rows.append(r)
M=sp.Matrix(rows); rk=M.rank(); print("rank of equations:",rk,"=> dim Pi_nt affine =",30-rk)
# one-sided 0/1 profiles in Pi_nt
cnt=0; classes={}
perms=list(itertools.permutations(range(5)))
def canon(F):
    best=None
    for p in perms:
        img=tuple(sorted(tuple(sorted(p[i] for i in B)) for B in F))
        if best is None or img<best: best=img
    return best
for choice in range(1<<15):
    F=[cuts[i] if (choice>>i)&1 else full-cuts[i] for i in range(15)]
    X=[sum(1 for B in F if a in B) for a in range(5)]
    if len(set(X))==1:
        cnt+=1; cl=canon(F); classes.setdefault(cl,0); classes[cl]+=1
print("one-sided near-tight profiles:",cnt,"classes:",len(classes))
# code vertices x_R=1[|B cap R| even], |R|=3, and complements
def xR(R): return [1 if len(B&R)%2==0 else 0 for B in Bs]
V=[]
for R in itertools.combinations(range(5),3):
    R=frozenset(R); v=xR(R); V.append(v); V.append([1-x for x in v])
for v in V:
    assert all(v[idx[B]]+v[idx[full-B]]==1 for B in Bs)
    X=[sum(v[idx[B]] for B in Bs if a in B) for a in range(5)]; assert len(set(X))==1
D=sp.Matrix([[v[i]-V[0][i] for i in range(30)] for v in V[1:]])
print("code hull affine dim:",D.rank())
avg=[sum(sp.Rational(v[i]) for v in V)/len(V) for i in range(30)]
print("uniform average == 1/2 everywhere:",all(x==sp.Rational(1,2) for x in avg))
x5=[1 if len(B)%2==0 else 0 for B in Bs]
D2=sp.Matrix([[v[i]-V[0][i] for i in range(30)] for v in V[1:]+[x5,[1-x for x in x5]]])
print("with even-set vertex and complement, affine dim:",D2.rank())
