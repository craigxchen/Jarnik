# Genuine parent with rank-0 tuples at inert l<=y: primes p=a^2+(Q b)^2, Q=prod_{l<=y, l=3(4)} l.
import math
from sympy import isprime, primerange
from tuples import gauss, omega, group_rank
y=23; ls=[l for l in primerange(3,y+1) if l%4==3]; Q=math.prod(ls)
pr=set()
for b in range(1,40):
    for a in range(1,400,2):
        p=a*a+(Q*b)**2
        if isprime(p): pr.add(p)
pr=sorted(pr)[:60]
vecs=[tuple(omega(gauss(p),l) for l in ls) for p in pr]
print("Q=",Q,"s=",len(pr),"min p=",pr[0],"Q^2=",Q*Q,"D=",len(set(vecs)),"rank=",group_rank(vecs,[l+1 for l in ls])[0])
print("min log p / theta3(y):", math.log(pr[0])/sum(math.log(l) for l in ls))
