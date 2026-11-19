from math import gcd,isqrt
def sq(n): return n>=0 and isqrt(n)**2==n
# homogenized: m=p/q. A*q^4 = -(p+3q)(p^3-11p^2q+39pq^2-13q^3); B*q^4=2(p^4+70p^2q^2-48pq^3+9q^4)
SA=[];SB=[];SmA=[];SmB=[]
for q in range(1,200):
    for p in range(-600,601):
        if gcd(p,q)!=1: continue
        a=-(p+3*q)*(p**3-11*p*p*q+39*p*q*q-13*q**3); b=2*(p**4+70*p*p*q*q-48*p*q**3+9*q**4)
        if a!=0 and sq(a): SA.append((p,q))
        if sq(b): SB.append((p,q))
        if a!=0 and sq(-a): SmA.append((p,q))
        if sq(-b): SmB.append((p,q))
print("A square:",SA[:15],len(SA)); print("B square:",SB[:15],len(SB)); print("-A sq:",SmA[:10],len(SmA)); print("-B sq:",SmB[:10],len(SmB))
