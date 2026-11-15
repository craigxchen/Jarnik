# Independent Gaussian-rational arithmetic and fractional gcd; C(t0) for the cube cluster.
from fractions import Fraction as Fr
from math import gcd, isqrt
import mpmath as mp
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def ggcd_int(a,b):
    # Gaussian integers as tuples of ints
    while b!=(0,0):
        n=b[0]*b[0]+b[1]*b[1]
        # a/b = a*conj(b)/n
        x=a[0]*b[0]+a[1]*b[1]; y=a[1]*b[0]-a[0]*b[1]
        q=((2*x+n)//(2*n),(2*y+n)//(2*n))
        qb=gmul(q,b); r=(a[0]-qb[0],a[1]-qb[1])
        a,b=b,r
    return a
def roots_cube(y,eps=(1,-1,-1,1)):
    r=[(Fr(e)/v, v) for e,v in zip(eps,y)]
    out=[]
    for x in r: out+=[x,(-x[0],-x[1])]
    return out
MASKS=[0b00000000,0b00001111,0b00110011,0b01010101,0b10101010,0b11001100,0b11110000,0b11111111]
def values(roots,t0,masks=MASKS):
    vals=[]
    for m in masks:
        v=(Fr(1),Fr(0))
        for k,(a,b) in enumerate(roots):
            if (m>>k)&1: b=-b
            v=gmul(v,(t0-a,-b))
        vals.append(v)
    return vals
def frac_gcd(vals):
    den=1
    for a,b in vals: den=den*a.denominator//gcd(den,a.denominator); den=den*b.denominator//gcd(den,b.denominator)
    ints=[(int(a*den),int(b*den)) for a,b in vals]
    g=(0,0)
    for z in ints: g=ggcd_int(z,g) if g!=(0,0) else z
    return g,den,ints   # gcd = g/den
def C_of(roots,t0,masks=MASKS,dps=40):
    vals=values(roots,t0,masks)
    g,den,ints=frac_gcd(vals)
    ng=g[0]**2+g[1]**2  # |g/den|^2 = ng/den^2
    # w_j = ints_j / g  exactly
    ws=[]
    for z in ints:
        x=z[0]*g[0]+z[1]*g[1]; yy=z[1]*g[0]-z[0]*g[1]
        assert x%ng==0 and yy%ng==0
        ws.append((x//ng,yy//ng))
    R2=ws[0][0]**2+ws[0][1]**2
    assert all(w[0]**2+w[1]**2==R2 for w in ws)
    assert len(set(ws))==len(ws)
    mp.mp.dps=max(dps,len(str(R2))+40)
    R=mp.sqrt(R2)
    # shortest arc containing all points: points close together; use max chord -> arc
    angs=sorted(mp.atan2(w[1],w[0]) for w in ws)
    gaps=[angs[i+1]-angs[i] for i in range(len(angs)-1)]+[angs[0]+2*mp.pi-angs[-1]]
    arc=(2*mp.pi-max(gaps))*R
    return arc/mp.sqrt(R), ws, R2, (g,den)
if __name__=="__main__":
    y=[Fr(-9,4),Fr(-2,3),Fr(1,4),Fr(8,3)]
    rts=roots_cube(y)
    for t0 in [Fr(13999,2),Fr(10**12+3999,2),Fr(10**30+423,2),Fr(10**30+424,2),Fr(10**30),Fr(10**30+1,3)]:
        C,ws,R2,_=C_of(rts,t0)
        print(t0 if t0.denominator==1 or len(str(t0))<30 else "big", mp.nstr(C,15), len(str(R2)))
