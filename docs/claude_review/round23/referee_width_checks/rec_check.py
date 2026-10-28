# Independent check of Lemmas 2.3-2.4 and the closed recursion (referee code; not reusing width_checks).
import sys
def H(n,x,y):
    lo,hi=min(x,y),max(x,y)
    if n==2: return 0
    if n==3: return lo
    if n==4: return min((x+y)//2, 2*lo)
    if n==5: return min(lo+hi//3, 3*lo)
def F(n,x,y): return x+y+H(n,x,y)
def peel(vals):
    v=sorted(vals,reverse=True); return max(v[i]+i for i in range(len(v)))
N=int(sys.argv[1]) if len(sys.argv)>1 else 60
bad=0
# monotonicity
for n in (2,3,4,5):
    for x in range(N+1):
        for y in range(N+1):
            if F(n,x+1,y)<F(n,x,y) or F(n,x,y+1)<F(n,x,y): bad+=1; print("mono fail",n,x,y)
# (R_n) and tightness
tight={3:0,4:0,5:0}
for n in (3,4,5):
    for a in range(N+1):
        for b in range(N+1):
            best=-1
            for p in range(a+1):
                for q in range(b+1):
                    X=H(n-1,a-p,b)+q; Y=H(n-1,a,b-q)+p
                    best=max(best,min(X,Y))
            if best>H(n,a,b): bad+=1; print("R fail",n,a,b)
            if best==H(n,a,b): tight[n]+=1
# direct: peel(V_n(a,b)) <= F_n(a,b) with F_{n-1} entries; and exact recursion rho
rho={2:{}}
M2=N
for a in range(M2+1):
    for b in range(M2+1): rho[2][(a,b)]=a+b
for n in (3,4,5):
    rho[n]={}
    for a in range(M2+1):
        for b in range(M2+1):
            V=[F(n-1,ap,b) for ap in range(a+1)]+[F(n-1,a,bp) for bp in range(b)]
            if peel(V)>F(n,a,b): bad+=1; print("peel fail",n,a,b)
            W=[rho[n-1][(ap,b)] for ap in range(a+1)]+[rho[n-1][(a,bp)] for bp in range(b)]
            rho[n][(a,b)]=peel(W)
            if rho[n][(a,b)]!=F(n,a,b): print("rho!=F",n,a,b,rho[n][(a,b)],F(n,a,b))
print("N=",N,"failures:",bad,"tight counts (H_n = max min):",tight,"of",(N+1)**2)
print("rho5(m,m) == floor(10m/3) for all m<=N:", all(rho[5][(m,m)]==10*m//3 for m in range(N+1)))
print("F5(m,m) == floor(10m/3):", all(F(5,m,m)==10*m//3 for m in range(10**4)))
print("F4(m,m)==3m, F3(m,m)==3m:", all(F(4,m,m)==3*m and F(3,m,m)==3*m for m in range(10**4)))
