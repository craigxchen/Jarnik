from check import *
def pmod(a,b):
    db=deg(b)
    while a and deg(a)>=db: a^=b<<(deg(a)-db)
    return a
def gcd2(a,b):
    while b: a,b=b,pmod(a,b)
    return a
def der(a): # derivative over F_2: keep odd-degree terms shifted
    r=0;i=1
    while (a>>i):
        if (a>>i)&1 and i%2==1: r|=1<<(i-1)
        i+=1
    return r
for r in (2,3):
    m=2**r;d=r+1;Q=2**d;L=(d-2)*2**(d-1);Y=((1<<L),1)
    N=(1,0)
    for v in range(2,Q): N=mul4(N,add4(Y,(v,0)))
    assert N[1]==0; N=N[0]
    print(m,"deg N",deg(N),"deg gcd(N,N')",deg(gcd2(N,der(N))))
