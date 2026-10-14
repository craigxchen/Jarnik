# Prop 5.4 check: Y_3(2), dim (J/grid)_7 vs L = grid + lines + quads + width powers of class-sum combos.
import itertools, sys, time
from fractions import Fraction
from rlib import *
from echelon import Ech
j=2; t=3
avec=[a for a in itertools.product((0,1),repeat=3) if any(a)]   # 7 coords
zs=list(itertools.product((0,1),repeat=3))
dot=lambda a,z: sum(x*y for x,y in zip(a,z))%2
def cz(y):  # class sums, times 4 (integer check)
    return [-sum(ya*(-1)**dot(a,z) for ya,a in zip(y,avec)) for z in zs]  # = 4*c_z/... c_z = -(1/4)sum y_a chi_a(z)*... check below
Y=[]
for y in itertools.product(range(-j,j+1),repeat=7):
    c4=cz(y)
    if all(v%4==0 for v in c4): Y.append(y)
print("|Y_3(2)| =",len(Y))
# sanity: recompute y from c
for y in Y[:50]:
    c=[v//4 for v in cz(y)]
    # c_z = -(1/4) sum_a y_a chi_a(z)  ->  check y_a = sum_{z: a.z=1} c_z  needs factor; fix scaling
    ys=[sum(cc for cc,z in zip(c,zs) if dot(a,z)) for a in avec]
    assert ys==[ -2*v for v in y] or ys==list(y), (y,ys)
print("scaling check", [sum(cc for cc,z in zip([v//4 for v in cz(Y[7])],zs) if dot(a,z)) for a in avec], Y[7])
cap=[2*j+1]*7
print("#monos deg<=7:",len(monos_le(7,7,cap)),"deg7:",len(monos_eq(7,7,cap)),"deg6:",len(monos_eq(7,6,cap)))
T0=time.time()
h={}
for e in (5,6,7):
    h[e]=h_func_exact(Y,e,cap); print("h_Y(%d) exact ="%e,h[e],time.time()-T0); sys.stdout.flush()
H7=h[7]-h[6]; H6=h[6]-h[5]
dJ7=len(monos_eq(7,7,cap))-H7; dJ6=len(monos_eq(7,6,cap))-H6
print("dim (J(Y)/grid)_6 =",dJ6," dim (J(Y)/grid)_7 =",dJ7)
import pickle; pickle.dump((Y,h,dJ6,dJ7),open('fano_Y.pkl','wb'))
