import math,random,numpy as np
exec(open('dirichlet.py').read().split('random.seed(1)')[0])
random.seed(1)
for Rp in [1e6,1e8]:
  for _ in range(5):
    th=mpf(random.random())*2*pi
    B=int(4*Rp**0.25)+5
    a=np.arange(-B,B+1); A,Bb=np.meshgrid(a,a); A=A.ravel(); Bb=Bb.ravel()
    m=(A*A+Bb*Bb>0)&(A*A+Bb*Bb<=B*B); A=A[m];Bb=Bb[m]
    nb=np.hypot(A,Bb); s=np.abs(np.sin(float(th)-np.arctan2(Bb,A))); v=nb*(1+math.sqrt(Rp)*s)
    i=np.argmin(v); print(Rp,float(th%(pi/2)), A[i],Bb[i],v[i], F_cf(th,Rp))
