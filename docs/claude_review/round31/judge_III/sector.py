# [J] check: #{w=a+bi: |w|^2<=X, |w|^2 odd >1, dist(arg w,(pi/4)Z)<=d} <= 8*tan(d)*(X+sqrt X) (no +O(sqrt X) boundary term)
import math,random
random.seed(1)
worst=0
for trial in range(60):
    X=random.choice([10**4,10**5,10**6]); d=random.choice([1e-1,1e-2,1e-3,3e-4])
    R=int(math.isqrt(X)); cnt=0
    for a in range(-R,R+1):
        for b in range(-R,R+1):
            n=a*a+b*b
            if n>X or n%2==0 or a*b==0: continue
            t=math.atan2(b,a); q=math.pi/4
            if abs(t-round(t/q)*q)<=d: cnt+=1
    bound=8*math.tan(d)*(X+math.sqrt(X))
    worst=max(worst,cnt/bound)
print('max ratio count/bound',worst)
