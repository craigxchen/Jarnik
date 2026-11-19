# Validate the (G1)->sigma reduction on an exact (degenerate) solution and on random NON-solutions.
import mpmath as mp, random
mp.mp.dps=40
E=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
inc=[[k for k,e in enumerate(E) if i in e] for i in range(5)]
def build(sig,c,K):
    a2=-c**3/K; q=c/a2; m=(1+q)/(1-q); L2=c/(m**2-1); mu=m*m
    def xi(i,s):
        v=s**-3*(s**2+a2)
        for k in inc[i]: v*=(s-sig[k])
        return v
    def parts(i,P):
        s=(-2*m*P+mp.sqrt(4*m*m*P*P+4*c))/2; si=-c/s; w=(c+s*s)/(2*s)
        return (xi(i,s)+xi(i,si))/2,(xi(i,s)-xi(i,si))/(2*w)
    return m,L2,mu,parts
def check(sig,c,K,label):
    m,L2,mu,parts=build(sig,c,K)
    Ps=[mp.mpf(3)+k for k in range(4)]
    out=[]
    for i in range(5):
        ys=[parts(i,P)[1] for P in Ps]
        xs=[(parts(i,P)[0]-P*parts(i,P)[1])/(P*P+L2) for P in Ps]
        out.append((mp.nstr(max(abs(y-ys[0]) for y in ys),3), mp.nstr(abs(xs[2]-2*xs[1]+xs[0]),3)))
    # n_i vanishes at rho_e for e at i
    van=[]
    for i in range(5):
        for k in inc[i]:
            s=sig[k]; P=(c-s*s)/(2*m*s)
            u,y=parts(i,P); x=(u-P*y)/(P*P+L2)
            X=x*P+y; Y=mp.sqrt(L2)*x
            van.append(abs(X*X+Y*Y-mu*y*y)/(abs(X*X)+1))
    print(label,' (y-spread, 2nd diff of x) per vertex:',out,' max |n_i(rho_e)|rel:',mp.nstr(max(van),3))
r=mp.mpf('0.2514'); c=mp.mpf(-1); K=r**4
sg=[1,1,-1,-1,-1,-1,1,1,-1,1]
check([s*r for s in sg],c,K,'exact degenerate solution')
random.seed(1)
check([mp.mpf(random.uniform(-2,2)) for _ in range(10)],c,mp.mpf(0.3),'random non-solution')
