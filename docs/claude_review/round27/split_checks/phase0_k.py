# Exact: k_j = s_j - s_min with s = d + n/d for the two five-divisor families; check z_j^2 - z_1^2 = 2 k_j s_1 + k_j^2.
from math import isqrt
def check(n,D,name):
    s=[d+n//d for d in D]; z=[d-n//d for d in D]
    i=min(range(5),key=lambda j:s[j]); s1,z1=s[i],z[i]
    k=[sj-s1 for sj in s]
    ok=all(z[j]**2-z1**2==2*k[j]*s1+k[j]**2 for j in range(5))
    print(name,"k =",sorted(k),"identity",ok)
u,w=1843,1064
for _ in range(4):
    n=u*(u+2)*(u+3)*(u+5)
    D=[(u+2)*(u+3),3*(w+1)*(w+3),(u+2)*(u+5),(u+3)*(u+5),3*(w+2)*(w+3)]
    check(n,D,"mnn31")
    u,w=2*u+3*w+7,u+2*w+4
X,Z=37,277
for _ in range(4):
    n=(X*X-1)*(X*X-4)*(X*X-16)*(X*X-49)
    F0=(X-1)*(X-2)*(X-4)*(X+7); F2=(X*X-4)*(X*X-16); F3=(X*X-1)*(X*X-16); F4=(X*X-1)*(X*X-4)
    F1=X**4-35*X*X+56+2*X*Z
    if n%F1: F1=X**4-35*X*X+56-2*X*Z
    D=[F0,F1,F2,F3,F4]
    assert all(n%d==0 for d in D)
    check(n,D,"m-czech")
    X,Z=15*X+2*Z,15*Z+112*X
