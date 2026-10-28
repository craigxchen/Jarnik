# Referee: constants of prove.md §§2-6, Lemma 3.1 brute force, Lemma 3.2 eigenvalue check.
import math, random
import numpy as np
pi=math.pi
print("Lemma3.1 consts (2+pi/2)/sqrt2=%.5f  pi/2=%.5f" % ((2+pi/2)/math.sqrt(2), pi/2))
print("Cor2.2 5*1.25506=%.5f" % (5*1.25506))
print("Cor4.3 32*sqrt(1/ln2)=%.4f 4*sqrt(1/ln2)=%.4f 40.4/pi=%.4f 39*sqrt(1.26)=%.4f 3*20.24=%.2f 32*1.26=%.3f" %
      (32*math.sqrt(1/math.log(2)),4*math.sqrt(1/math.log(2)),40.4/pi,39*math.sqrt(1.26),3*20.24,32*1.26))
print("Prop6.1 (4/pi)(1+2^-1/2)=%.5f  atan(1/2)+atan(1/3)-pi/4=%.2e  asin(5^-1/2)-atan(1/2)=%.2e" %
      ((4/pi)*(1+2**-0.5), math.atan(.5)+math.atan(1/3)-pi/4, math.asin(5**-0.5)-math.atan(.5)))
print("16 sqrt2=%.4f" % (16*math.sqrt(2)))
x=math.exp(-5); print("tail ratio at w=10: %.5f" % ((4/pi)*(math.asin(x)+math.asin(x/math.sqrt(2)))/x))

# Lemma 3.1 brute force: associate classes alpha != 0, Norm <= Y, arg mod pi/2 in arc I of length eta
random.seed(1)
def classes(Y):
    out=[]
    r=int(math.isqrt(Y))
    for a in range(1,r+1):          # representatives with a>0, b>=0 (first quadrant incl. positive real axis)
        for b in range(0,r+1):
            if a*a+b*b<=Y: out.append(math.atan2(b,a))
    return np.array(sorted(out))
worst=0
for Y in [50,200,1000,5000,20000]:
    ang=classes(Y); m=len(ang)
    for trial in range(400):
        eta=random.choice([random.uniform(0,pi/2), random.uniform(0,10/Y), random.uniform(0,1/math.sqrt(Y))])
        c=random.uniform(0,pi/2)
        # count in arc [c, c+eta] mod pi/2 (closed)
        t=(ang-c)%(pi/2)
        cnt=int(np.sum(t<=eta+1e-15))
        bd=eta*Y/2+(2+eta)/math.sqrt(2)*math.sqrt(Y)+pi/2
        worst=max(worst,cnt/bd)
    # also worst arcs anchored at a point
    for eta in [1/Y, 3/Y, 1/math.sqrt(Y)]:
        for i in range(0,m,max(1,m//300)):
            t=(ang-ang[i])%(pi/2); cnt=int(np.sum(t<=eta+1e-15))
            bd=eta*Y/2+(2+eta)/math.sqrt(2)*math.sqrt(Y)+pi/2; worst=max(worst,cnt/bd)
print("Lemma 3.1 worst count/bound ratio = %.4f (must be <=1)" % worst)

# Lemma 3.2: optimal constant Lambda = top eigenvalue of G_{nm} = sum_{|k|<=K} e(k(x_n-x_m))
def Lam(x,K):
    d=x[:,None]-x[None,:]
    G=np.ones_like(d)
    for k in range(1,K+1): G=G+2*np.cos(2*pi*k*d)
    return np.linalg.eigvalsh(G)[-1]
def maxarc(x,length,closed=False):
    xs=np.sort(x%1.0); n=len(xs); ext=np.concatenate([xs,xs+1]); best=0
    for i in range(n):
        j=np.searchsorted(ext, xs[i]+length, side='right' if closed else 'left')
        best=max(best,j-i)
    return best
bad=0
for trial in range(60):
    n=random.choice([20,60,150]); K=random.choice([2,5,12,30])
    mode=trial%3
    if mode==0: x=np.random.rand(n)
    elif mode==1: x=np.concatenate([np.random.rand(n//2), 0.3+0.01*np.random.rand(n-n//2)])
    else: x=(np.arange(n)/n + 1e-3*np.random.rand(n))
    L=Lam(x,K)
    delta=1/(2*K); m=maxarc(x,delta)
    up=2*m*(2*K+1/delta)
    t=maxarc(x,1/(8*K),closed=True); lo=(2*K+1)*t/math.sqrt(2)
    if not (lo<=L*(1+1e-9) and L<=up*(1+1e-9)): bad+=1; print("VIOL",n,K,mode,lo,L,up)
print("Lemma 3.2 two-sided check: violations =",bad)
