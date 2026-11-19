# Independent sympy check of one F_13 toy example with 5 isotropic blocks and deg delta = 2.
import sympy as sp, sys
from itertools import combinations
T=sp.symbols('T'); p=13
line=[l for l in open('ex13.txt') if 'niso=5' in l][0]
cols=[tuple(map(int,s.strip('()').replace(';',',').split(','))) for s in line.split('|')[1].split()]
F=lambda c: sp.Poly(c,T,modulus=p)
V=[(F(T**2+b*T+g),F(T**2+bp*T+gp)) for b,g,bp,gp in cols]
roots={}
for i,j in combinations(range(5),2):
    d=V[i][0]*V[j][1]-V[i][1]*V[j][0]; assert d.degree()==1, d
    roots[(i,j)]=(-d.all_coeffs()[1]*pow(int(d.all_coeffs()[0]),-1,p))%p
assert len(set(roots.values()))==10
for i in range(5):
    assert sp.gcd(V[i][0],V[i][1]).degree()==0
    n=F(1)
    for e,a in roots.items():
        if i in e: n=n*F(T-a)
    dl=V[i][0]*V[i][1]-n
    iso=[(e[0]+1,e[1]+1) for e,a in roots.items() if i in e and (V[i][0]*V[i][1]).eval(a)%p==0]
    print('mark',i+1,'deg delta',dl.degree(),'isotropic blocks',iso)
