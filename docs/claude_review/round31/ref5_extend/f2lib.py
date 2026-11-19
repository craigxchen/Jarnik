# Own F_2[t] (ints as bit vectors) and F_4[t] = F_2[t][theta]/(theta^2+theta+1) arithmetic.
def deg(a): return a.bit_length()-1  # deg 0 = -1
def mul(a,b):
    if a.bit_length()<b.bit_length(): a,b=b,a
    r=0
    while b:
        lb=b&-b; r^=a<<(lb.bit_length()-1); b^=lb
    return r
def divmod2(a,b):
    q=0; db=deg(b)
    while deg(a)>=db:
        s=deg(a)-db; q|=1<<s; a^=b<<s
    return q,a
def gcd(a,b):
    while b: a,b=b,divmod2(a,b)[1]
    return a
def sq(a):
    r=0;i=0
    while a:
        if a&1: r|=1<<(2*i)
        a>>=1;i+=1
    return r
# F_4[t]: (A,B) = A + B*theta
def m4(x,y):
    A,B=x;C,D=y
    ac=mul(A,C); bd=mul(B,D)
    return (ac^bd, mul(A,D)^mul(B,C)^bd)
def add4(x,y): return (x[0]^y[0],x[1]^y[1])
def conj4(x): return (x[0]^x[1],x[1])
def norm4(x):
    A,B=x; return sq(A)^mul(A,B)^sq(B)
def deg4(x): return max(deg(x[0]),deg(x[1]))
def prod4(lst):
    r=(1,0)
    for x in lst: r=m4(r,x)
    return r
def lc4(x):
    d=deg4(x); return ((x[0]>>d)&1,(x[1]>>d)&1)
INV4={(1,0):(1,0),(0,1):(1,1),(1,1):(0,1)}
def shift4(x,s): return (x[0]<<s,x[1]<<s)
def mod4(x,y):
    dy=deg4(y); iy=INV4[lc4(y)]
    while deg4(x)>=dy and (x[0] or x[1]):
        c=m4(lc4(x),iy); x=add4(x,shift4(m4(c,y),deg4(x)-dy))
    return x
def gcd4(x,y):
    while y[0] or y[1]: x,y=y,mod4(x,y)
    if x[0] or x[1]: x=m4(INV4[lc4(x)],x)
    return x
