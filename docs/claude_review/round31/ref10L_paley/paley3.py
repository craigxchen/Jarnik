# sizes: full-core condition, 5/8 exponent, deletion threshold, with real consecutive primes from P
import math
def isp(n): return n>1 and all(n%d for d in range(2,int(n**.5)+1))
def primes_from(P,cnt):
    out=[];n=P
    while len(out)<cnt:
        if isp(n): out.append(n)
        n+=1
    return out
def mlist(lo,hi): return [m for m in range(lo,hi) if isp(m) and m%4==3]
for n in (1,2):
  for C in (1,10):
    for m in mlist(7,400):
        P=max(int(m*math.log(m)**2),3)
        pr=primes_from(P,n*m)
        while pr[-1]>2*P: P=int(P*1.1)+1; pr=primes_from(P,n*m)
        lq=sum(math.log(p) for p in pr)
        Q=[sum(math.log(p) for p in pr[n*v:n*v+n]) for v in range(m)]
        Qs=sorted(Q)
        lg_min=sum(Qs[:(m-3)//4])  # lower bound for min triple g (actual min over triples >= this)
        full = (lg_min>2*math.log(C)) and (math.log((m+1)*C/2)+Qs[-1]+Qs[-2] < lq/2)
        if full:
            print("n=%d C=%d first m with full-core exclusion (P=m log^2 m=%d): m=%d"%(n,C,P,m)); break
for n in (1,2):
    for m in [43,103,503,1019,5003,20011]:
        if not isp(m): continue
        P=int(m*math.log(m)**2); pr=primes_from(P,n*m); lq=sum(math.log(p) for p in pr)
        lQ=n*math.log(2*P)
        # her inequality: (r+1) n log(2P) + r/2 log r + log(C m (r+1)) < lq/2, C=1, r=k+1
        kmax=-1
        for k in range(m):
            r=k+1
            if (r+1)*lQ+0.5*r*math.log(r)+math.log(m*(r+1)) < lq/2: kmax=k
            else: break
        th=math.log(m)/math.log(P)
        print("n=%d m=%d P=%d: kmax/m=%.3f ; predicted n/(2n+theta)=%.3f (theta=log m/log P=%.3f); limit n/(2n+1)=%.3f"%(n,m,P,kmax/m,n/(2*n+th),th,n/(2*n+1)))
# 5/8: log g/log q for actual primes (equal n per column)
for n in (1,2):
    for m in [43,503,5003]:
        P=int(m*math.log(m)**2); pr=primes_from(P,n*m); lq=sum(math.log(p) for p in pr)
        Q=sorted(sum(math.log(p) for p in pr[n*v:n*v+n]) for v in range(m))
        print("n=%d m=%d: min log g/log q >= %.4f ; 1/2+that/2 = %.4f"%(n,m,sum(Q[:(m-3)//4])/lq,0.5+sum(Q[:(m-3)//4])/lq/2))
