import numpy as np, random, math
rng = random.Random(5)
def walsh(v):
    v=v.astype(np.float64).copy();h=1
    while h<len(v):
        for i in range(0,len(v),2*h):
            a=v[i:i+h].copy();b=v[i+h:i+2*h].copy();v[i:i+h]=a+b;v[i+h:i+2*h]=a-b
        h*=2
    return v
for p in (10007,100003,1000003):
    u=np.arange(p,dtype=np.int64); chi=-np.ones(p,dtype=np.int8); chi[(u*u)%p]=1
    S=[];DEV=[];byk={2:[],4:[],6:[],8:[]}
    for _ in range(300):
        xs=[rng.randrange(p) for _ in range(8)]
        R=np.stack([chi[(u+x)%p] for x in xs])
        code=((R[1:]!=R[0]).astype(np.int64)*(1<<np.arange(7))[:,None]).sum(0)
        pat=np.bincount(code,minlength=128)
        m=walsh(pat)/p
        for A in range(1,128):
            k=bin(A).count('1'); k=k+1 if k%2 else k
            byk[k].append(m[A]**2)
        S.append(np.abs(m[1:]).sum()); DEV.append(np.abs(pat*128-p).max()/p)
    print(p,'mean sum|m_A|=%.4f  *sqrt(p)/127=%.3f  mean maxdev=%.4f'%(np.mean(S),np.mean(S)*math.sqrt(p)/127,np.mean(DEV)),
          ' p*E m_A^2 by |A|:',{k:round(np.mean(v)*p,3) for k,v in byk.items()})
