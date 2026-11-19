import numpy as np
from common import *
rng=np.random.default_rng(7)
for tr in range(4):
    th=rng.uniform(-np.pi,np.pi,5); th-=th[4]
    t,n,P,_=real_sample(rng,theta=th)
    p=P; a=np.conj(P); q=p/a   # q5=1 since theta5=0
    d={}
    for (i,j) in E: d[(i,j)]=(q[i]-q[j])/(-2j*t[eid(i,j)]); d[(j,i)]=(q[j]-q[i])/(-2j*(-t[eid(i,j)]))
    K=np.array([np.prod([d[(i,j)] for j in range(5) if j!=i]) for i in range(5)])
    B=np.prod(a)**2
    e1=abs(B**3-np.prod(q)/np.prod(K))/abs(B**3)
    e2=max(abs(n[eid(i,j)]-d[(i,j)]*q[i]*q[j]/(K[i]*K[j]*B))/n[eid(i,j)] for (i,j) in E)
    print("rel err B^3 eq",e1," n formula",e2)
