# Independent symbolic checks for the four-block star (Part V.3)
import sympy as sp
from itertools import combinations
v,t,k=sp.symbols('v t k')
r=t**2; f=(1-r)*(r**2-3*r+1)
x=(k+f/k)/2; y=(f/k-k)/(2*t)
w=[v*(1+x)/2, v*(1-x)/2, v*(-r+y)/2, v*(-r-y)/2]
e=[sp.together(sum(sp.prod(c) for c in combinations(w,m))) for m in (1,2,3,4)]
print('s3-s1s2 identically zero:', sp.simplify(e[2]-e[0]*e[1])==0)
# generic G,A identity with symbolic w satisfying s3=s1 s2 : check with symbols
W=sp.symbols('w2:6'); T=sp.symbols('T')
S=[sum(sp.prod(c) for c in combinations(W,m)) for m in (1,2,3,4)]
G=sp.expand(sp.prod([T-sp.I*wj for wj in W]))
A=-T/S[0]+sp.I
br=lambda F,H: sp.expand(sp.re(sp.expand(F))*sp.im(sp.expand(H))-sp.im(sp.expand(F))*sp.re(sp.expand(H)))
TT=sp.symbols('TT',real=True)
Wr=sp.symbols('a2:6',real=True)
sub={W[m]:Wr[m] for m in range(4)}
Gr=G.subs(sub).subs(T,TT); Ar=A.subs(sub).subs(T,TT); Sr=[s.subs(sub) for s in S]
d=sp.simplify(br(Gr,Ar)-Sr[3]-(Sr[2]/Sr[0]-Sr[1])*TT**2)
print('[G,A]=s4+(s3/s1-s2)T^2 :', d==0)
# Jacobian rank of direction map at a second, independent point
D=[wj**-2-1 for wj in w]
for pt in [{v:sp.Rational(3,2),t:2,k:2},{v:sp.Rational(5,7),t:3,k:sp.Rational(-4,3)},{v:2,t:sp.Rational(1,3),k:5}]:
    M=sp.Matrix([[Dj.subs(pt)]+[sp.diff(Dj,z).subs(pt) for z in (v,t,k)] for Dj in D])
    ws=[wj.subs(pt) for wj in w]
    print(pt,'w=',ws,'det=',M.det())
