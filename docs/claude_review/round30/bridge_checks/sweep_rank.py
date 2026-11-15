# Dimension of the subvariety of T_5 swept by the 5-subset curves of the cube family:
# Jacobian rank of (y0,y1,t) -> log-ratio coordinates (z_j/z_0), j=1..4, at random real points (complex step-free, mpmath).
import itertools, mpmath as mp
mp.mp.dps=40
names=['0','1','A','Ac','B','Bc','P','Pc']
sets={'0':[], '1':list(range(8)), 'A':[0,1,2,3],'Ac':[4,5,6,7],'B':[0,1,4,5],'Bc':[2,3,6,7],'P':[0,2,4,6],'Pc':[1,3,5,7]}
def ycomplete(y0,y1):
    # y2 from y0 y1 y2 (y0+y1+y2) = -1  : y0y1 y2^2 + y0y1(y0+y1) y2 + 1 = 0
    A=y0*y1; B=A*(y0+y1); disc=B*B-4*A
    y2=(-B+mp.sqrt(disc))/(2*A); return [y0,y1,y2,-(y0+y1+y2)]
def zvals(y0,y1,t):
    y=ycomplete(y0,y1); eps=[1,-1,-1,1]
    r=[e/v+1j*v for e,v in zip(eps,y)]
    roots=[]
    for x in r: roots+=[x,-x]
    out={}
    for nm,F in sets.items():
        p=mp.mpc(1)
        for k in range(8):
            rho=mp.conj(roots[k]) if k in F else roots[k]
            p*= (t-rho)
        out[nm]=p
    return out
def coords(S,v):
    y0,y1,t=v; z=zvals(y0,y1,t)
    c=[]
    for s in S[1:]:
        w=mp.log(z[s]/z[S[0]])   # norm-one ratio: log is i*angle; keep real-imag
        c+= [mp.re(w), mp.im(w)]
    return c
if __name__=="__main__":
    import random
    random.seed(3)
    from collections import Counter
    out=Counter()
    for S in itertools.combinations(names,5):
        ranks=[]
        for trial in range(2):
            while True:
                y0=mp.mpf(random.uniform(-3,-0.5)); y1=mp.mpf(random.uniform(-1,-0.2))
                A=y0*y1; B=A*(y0+y1)
                if B*B-4*A>0: break
            v=[y0,y1,mp.mpf(random.uniform(1,3))]
            h=mp.mpf(10)**-12
            J=[]
            for i in range(3):
                vp=list(v); vm=list(v); vp[i]+=h; vm[i]-=h
                J.append([(a-b)/(2*h) for a,b in zip(coords(S,vp),coords(S,vm))])
            M=mp.matrix(J)
            sv=mp.svd_r(M,compute_uv=False)
            ranks.append(sum(1 for s in sv if s>mp.mpf(10)**-15))
        out[(tuple(ranks))]+=1
    print("Jacobian ranks of (y0,y1,t) -> T_5 over the 56 five-subsets:",out)
