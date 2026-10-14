import time, pickle, itertools, math
import sympy as sp
from m8lib import *
A, cols, H, JA, W = pickle.load(open('m8_step1.pkl','rb'))
gw = pickle.load(open('m8_gw.pkl','rb'))
a,b,c = 1,2,3
ca, cb, cc = cols[a], cols[b], cols[c]
z = tuple((x+y+w)//2 for x,y,w in zip(ca,cb,cc))
print("c_a,c_b,c_{a+b} =", ca, cb, cc, " z = (c_a+c_b+c_c)/2 =", z, " parity ok:", all((x+y+w)%2==0 for x,y,w in zip(ca,cb,cc)))
for F, name in [((ca,cb,cc),'(y_a,y_b,y_{a+b})'), ((ca,cb,z),'(y_a,y_b,z)')]:
    g, rho, ny = proj_ideal_gens(A, F, 8)
    print(name, "|Y|", ny, "rho", rho)
    Y = sp.symbols('y1 y2 y3')
    for e, fl in g.items():
        for f in fl:
            expr = sum(cf*sp.prod([Y[i]**ex[i] for i in range(3)]) for ex, cf in f.items())
            print("  deg", e, ":", sp.factor(expr))
