import random
from f24 import *
random.seed(1)
def rnd4(d):
    return (random.getrandbits(d+1), random.getrandbits(d+1))
for _ in range(200):
    P=rnd4(random.randint(0,30)); Q=rnd4(random.randint(0,30)); R=rnd4(random.randint(0,30))
    assert mul4(P,add4(Q,R))==add4(mul4(P,Q),mul4(P,R))
    assert mul4(mul4(P,Q),R)==mul4(P,mul4(Q,R))
    assert norm4(P)==mul4(P,conj4(P))[0] and mul4(P,conj4(P))[1]==0
    if not iszero4(Q):
        q,r=divmod4(P,Q); assert add4(mul4(q,Q),r)==P and (iszero4(r) or deg4(r)<deg4(Q))
    g=gcd4(mul4(P,R),mul4(Q,R))
    if not iszero4(R): 
        _,r=divmod4(g,monic4(R)); assert iszero4(r)
    a=random.getrandbits(40); b=random.getrandbits(40)
    q,r=divmod2(a,b) if b else (0,a)
    if b: assert mul2(q,b)^r==a
# theta^2=theta+1
th=(0,1); assert mul4(th,th)==(1,1)
print("ok")
