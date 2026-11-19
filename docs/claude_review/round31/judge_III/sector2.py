import math
def count(X,d):
    R=math.isqrt(X); c=0
    for a in range(-R,R+1):
        for b in range(-R,R+1):
            n=a*a+b*b
            if n>X or n%2==0 or a*b==0: continue
            t=math.atan2(b,a); q=math.pi/4
            if abs(t-round(t/q)*q)<=d: c+=1
    return c
for X in [10**4,10**5]:
    for d in [1e-1,1e-2,1e-3,3e-4]:
        c=count(X,d); print(X,d,c, round(c/(math.tan(d)*X),2))
