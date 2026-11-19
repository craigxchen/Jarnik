from mpmath import mp, atan2, mpf
from math import prod
from itertools import combinations
exec(open('inst.py').read().split("s=-K")[0].replace("print(","(lambda *a,**k:None)("))
mp.dps=80
ang={i:2*atan2(V[i][1],V[i][0]) for i in V}
diam=max(ang.values())-min(ang.values())
print("diam*N^(3/10)=",diam*mpf(N)**mpf(0.3), "2*max|R_i-R_j|=",2*(max(R.values())-min(R.values())))
