# Krawczyk certificate for a 5-point equidistant cluster at D=8 (class from equi5_found.pkl).
import pickle, numpy as np, mpmath as mp, itertools
from fractions import Fraction
ci,dl,pts,x=pickle.load(open('equi5_found.pkl','rb'))
D=8;L=3
mp.mp.dps=60; iv=mp.iv; iv.dps=60
# fixed: x0=0,y0=1, x1,y1 rounded to rationals; unknowns: x2..x7,y2..y7 (12)
fx1=Fraction(round(x[1]*1000),1000); fy1=Fraction(round(x[D+1]*1000),1000)
def build(u,ctx):
    X=[ctx.mpf(0),ctx.mpf(fx1.numerator)/fx1.denominator]+list(u[:6])
    Y=[ctx.mpf(1),ctx.mpf(fy1.numerator)/fy1.denominator]+list(u[6:])
    return X,Y
def imp(xk,yk,m):
    # Im (x+iy)^m expanded as real polynomial
    return [yk, 2*xk*yk, 3*xk**2*yk-yk**3][m-1]
def Fv(u,ctx):
    X,Y=build(u,ctx)
    return [sum(int(dl[i][k])*imp(X[k],Y[k],m) for k in range(D) if dl[i][k]) for i in range(4) for m in (1,2,3)]
def J(u,ctx):
    X,Y=build(u,ctx); rows=[]
    for i in range(4):
        for m in (1,2,3):
            r=[]
            for k in range(2,D):
                dx=[0,2*Y[k],6*X[k]*Y[k]][m-1]
                r.append(int(dl[i][k])*dx)
            for k in range(2,D):
                dy=[1,2*X[k],3*X[k]**2-3*Y[k]**2][m-1]
                r.append(int(dl[i][k])*dy)
            rows.append(r)
    return rows
# Newton refine
u=[mp.mpf(v) for v in list(x[2:D])+list(x[D+2:])]
# first a few damped LS steps to land on slice
for it in range(60):
    Fn=mp.matrix(Fv(u,mp)); Jn=mp.matrix(J(u,mp))
    du=mp.lu_solve(Jn,Fn); u=[a-b for a,b in zip(u,du)]
print('residual after Newton',mp.norm(mp.matrix(Fv(u,mp))))
Y_=mp.inverse(mp.matrix(J(u,mp)))
r=mp.mpf(10)**-30
Xb=[iv.mpf([a-r,a+r]) for a in u]
xt=[iv.mpf(a) for a in u]
Fx=Fv(xt,iv); JX=J(Xb,iv)
n=12
K=[]
for i in range(n):
    s=xt[i]-sum(iv.mpf(Y_[i,j])*Fx[j] for j in range(n))
    for j in range(n):
        cij=(1 if i==j else 0)-sum(iv.mpf(Y_[i,l])*JX[l][j] for l in range(n))
        s+=cij*(Xb[j]-xt[j])
    K.append(s)
inside=all((K[i].a>Xb[i].a) and (K[i].b<Xb[i].b) for i in range(n))
print('Krawczyk K(X) in interior of X:',inside, ' max width K',max(float(k.delta) for k in K))
# nondegeneracy on the box
Xr,Yr=build(Xb,iv)
roots=[(Xr[k],Yr[k]) for k in range(D)]
allr=roots+[(a,-b) for a,b in roots]
sep=min(float((abs(a1-a2)+abs(b1-b2)).a) for (a1,b1),(a2,b2) in itertools.combinations(allr,2) )
print('min |Im rho| lower bound',min(float(abs(b).a) for a,b in roots),' min L1 separation of the 16 roots (lower est)',sep)
for k in range(D): print('rho_%d ='%k, mp.nstr(build(u,mp)[0][k],12), '+ i', mp.nstr(build(u,mp)[1][k],12))
print('class',ci,'points (s1=1 frame):'); print(pts)
# full fibre profile on the certified box: points s in {+-1}^8 relative to the s1=1 frame (rho_k as computed, unnormalised)
mom=[[imp(Xr[k],Yr[k],m) for m in (1,2,3)] for k in range(D)]
allp=list(itertools.product([1,-1],repeat=D))
mv=[[sum(s_[k]*mom[k][m] for k in range(D)) for m in range(3)] for s_ in allp]
base=allp.index(tuple([1]*D))
same=[];undec=0
for i,v in enumerate(mv):
    z=[v[m]-mv[base][m] for m in range(3)]
    if all(zz.a<=0<=zz.b for zz in z): same.append(i)
print('points whose moment-difference box contains 0 (fibre of s1, upper bound):',len(same))
# global profile upper bound: group by overlapping boxes
import collections
cnt=collections.Counter()
for i,v in enumerate(mv):
    c=0
    for j,w in enumerate(mv):
        if all((v[m]-w[m]).a<=0<=(v[m]-w[m]).b for m in range(3)): c+=1
    cnt[c]+=1
print('max over all base points of #points not certified-distinct:',max(cnt))
pickle.dump((ci,dl,pts,[str(a) for a in u],fx1,fy1),open('equi5_cert.pkl','wb'))
