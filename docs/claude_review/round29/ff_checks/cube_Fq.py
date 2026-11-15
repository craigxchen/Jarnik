# Cube 8-point clusters over F_q(t), q = 3 mod 4, p>=5: z in F_{q^2}[t], conj = Frobenius (i -> -i).
import itertools, sys
def run(q,maxtries=50):
    def mul(a,b): return ((a[0]*b[0]-a[1]*b[1])%q,(a[0]*b[1]+a[1]*b[0])%q)
    def sub(a,b): return ((a[0]-b[0])%q,(a[1]-b[1])%q)
    def neg(a): return ((-a[0])%q,(-a[1])%q)
    def cj(a): return (a[0],(-a[1])%q)
    def pmul(P,Qp):
        R=[(0,0)]*(len(P)+len(Qp)-1)
        for i,a in enumerate(P):
            for j,b in enumerate(Qp):
                m=mul(a,b); R[i+j]=((R[i+j][0]+m[0])%q,(R[i+j][1]+m[1])%q)
        return R
    def lin(r): return [neg(r),(1,0)]   # t - r, low degree first
    inv=lambda a: pow(a,q-2,q)
    found=[]
    for a,b,c in itertools.product(range(1,q),repeat=3):
        d=(-(a+b+c))%q
        if d==0 or (a*b*c*d)%q!=1: continue
        y=[a,b,c,d]
        if any((y[i]+y[j])%q==0 for i in range(4) for j in range(i+1,4)): continue
        eps=[1,-1,-1,1]
        r=[((e*inv(v))%q, v) for e,v in zip(eps,y)]
        roots=[r[0],neg(r[0]),r[1],neg(r[1]),r[2],neg(r[2]),r[3],neg(r[3])]
        allr=roots+[cj(x) for x in roots]
        if len(set(allr))<16: continue
        Q=[(1,0)]
        for x in roots: Q=pmul(Q,lin(x))
        cnt=0
        for mask in range(256):
            z=[(1,0)]
            for k,x in enumerate(roots): z=pmul(z,lin(cj(x) if (mask>>k)&1 else x))
            diff=[sub(u,v) for u,v in zip(z,Q)]
            deg=max([i for i,u in enumerate(diff) if u!=(0,0)],default=-1)
            if deg<=4: cnt+=1
        found.append((tuple(y),cnt))
        if len(found)>=maxtries: break
    return found
for q in (7,11,19,23,31,43):
    f=run(q,8)
    print(q, "nondegenerate y's tested:",len(f)," cluster sizes:",sorted(set(c for _,c in f)), f[:2])
