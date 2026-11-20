# [J] Test 1: inert prime-power collision bound proposed by stepanov/cm/lifts analysts:
#   sum_{i<j} v_l(z_i - z_j) >= sum_{k>=1} max(0, M^2/(2(l+1)l^(k-1)) - M/2)
# for M distinct lattice points of |z|^2=N, l = 3 mod 4, l not dividing N.
# Also: class count of {z mod l^k : z zbar = N mod l^k} equals (l+1) l^(k-1).
import itertools, random, math
from math import log

def pts(N):
    out=[]; a=0
    while a*a<=N:
        b2=N-a*a; b=math.isqrt(b2)
        if b*b==b2:
            for x,y in {(a,b),(-a,b),(a,-b),(-a,-b),(b,a),(-b,a),(b,-a),(-b,-a)}:
                out.append((x,y))
        a+=1
    return sorted(set(out), key=lambda p: math.atan2(p[1],p[0]))

def vl(x,y,l):  # l-adic valuation of x+iy for inert l = valuation of gcd(x,y)
    if x==0 and y==0: return 10**9
    v=0
    while x%l==0 and y%l==0: x//=l; y//=l; v+=1
    return v

def classcount(N,l,k):
    m=l**k
    return sum(1 for x in range(m) for y in range(m) if (x*x+y*y-N)%m==0)

for l in (3,7,11):
    for k in (1,2):
        for N in (5*13*17, 65*29):
            c=classcount(N,l,k); assert c==(l+1)*l**(k-1),(l,k,N,c)
print("class counts (l+1)l^(k-1): OK for l in 3,7,11, k in 1,2")

def bound(M,l):
    s=0.0;k=1
    while True:
        t=M*M/(2*(l+1)*l**(k-1))-M/2
        if t<=0: break
        s+=t;k+=1
    return s

random.seed(1)
Ns=[5*13*17*29*37, 5*13*17*29*37*41, 5**2*13*17*29*41, 65*85*29*37*53]
inert=[3,7,11,19,23,31,43,47]
viol=0; tests=0; worst=1e9
for N in Ns:
    P=pts(N)
    for trial in range(300):
        M=random.choice([6,10,16,24,32])
        if M>len(P): continue
        if trial%2==0:
            st=random.randrange(len(P)); S=[P[(st+i)%len(P)] for i in range(M)]  # consecutive in angle
        else:
            S=random.sample(P,M)
        for l in inert:
            if N%l==0: continue
            lhs=sum(vl(a[0]-b[0],a[1]-b[1],l) for a,b in itertools.combinations(S,2))
            rhs=bound(M,l); tests+=1
            if lhs<rhs-1e-9: viol+=1
            worst=min(worst,lhs-rhs)
print("tests",tests,"violations",viol,"min slack",round(worst,3))

# level>=2 contribution constant: sum_{l=3 mod 4} log l/(l^2-1)
def primes(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
c=sum(log(p)/(p*p-1) for p in primes(10**6) if p%4==3)
print("sum_{l=3(4)} log l/(l^2-1) =",round(c,4))
# effect on W: Theorem 1.1 is (M^2/2)*(stuff)*log; extra M^2/2*c in sum log|c_ij|; via cut identity
# (8 * sum log|c_ij|) on the 8th-power side vs W*M^2-ish => O(M) in W. Report ratio extra/(level-1 mass) at M:
for M in (10,100,1000):
    lvl1=sum(log(p)/(p+1) for p in primes(M) if p%4==3)
    print("M",M,"level-1 inert mass",round(lvl1,3),"levels>=2 add",round(c,3))
