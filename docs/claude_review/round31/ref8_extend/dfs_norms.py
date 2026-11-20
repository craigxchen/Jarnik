import sys
SPLIT=[p for p in range(5,1000) if p%4==1 and all(p%d for d in range(2,int(p**.5)+1))]
def gen(NMAX,minpts,PMAX):
    out=[]
    P=[p for p in SPLIT if p<=PMAX]
    def rec(i,N,fac,t):
        if t>=minpts: out.append(dict(fac))
        for j in range(i,len(P)):
            p=P[j]
            if N*p>NMAX: break
            e=1; M=N*p
            while M<=NMAX:
                fac[p]=e; rec(j+1,M,fac,t*(e+1)); del fac[p]
                e+=1; M*=p
    rec(0,1,{},1)
    return out
if __name__=="__main__":
    L=gen(10**int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]))
    print(len(L))
