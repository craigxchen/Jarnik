# (i) block-norm gcds in mixed-rate family over many n and seeds; (ii) Gaussian gcd of the 4 points;
# (iii) bound [v2,v3][v2,Dv3] for odd part of gcd(N2,N3)
import random, math
from pell4 import *
random.seed(7)
def cross(a,b): return a[0]*b[1]-a[1]*b[0]
def D(v): return (v[0]+v[1], v[0]-v[1])
def odd(m):
    while m%2==0: m//=2
    return m
worst={}
for trial in range(40):
    while True:
        seeds=[(random.randint(-30,30),random.randint(-30,30)) for _ in range(3)]
        if all(math.gcd(*v)==1 for v in seeds) and cross(seeds[1],seeds[2])*cross(seeds[1],D(seeds[2]))!=0: break
    c=abs(cross(seeds[1],seeds[2])*cross(seeds[1],D(seeds[2])))
    mx12=mx13=mx23=mxg=1; viol=0
    for n in range(1,41):
        bl=[Apow(v,r*n) for v,r in zip(seeds,(6,7,7))]
        Ns=[N(b) for b in bl]
        g12,g13,g23=math.gcd(Ns[0],Ns[1]),math.gcd(Ns[0],Ns[2]),math.gcd(Ns[1],Ns[2])
        if c % odd(g23): viol+=1
        pts=words(bl); g=pts[0]
        for p in pts[1:]: g=ggcd(g,p)
        mx12=max(mx12,g12); mx13=max(mx13,g13); mx23=max(mx23,g23); mxg=max(mxg,N(g))
    print(seeds, "max gcd(N1,N2)",mx12,"max gcd(N1,N3)",mx13,"max gcd(N2,N3)",mx23,"cross bound",c,"violations",viol,"max N(gausgcd)",mxg)
