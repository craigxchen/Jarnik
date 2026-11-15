import sys, numpy as np, glob
sys.path.insert(0,'.'); from certify import certify; from basesolve import degen
from fullcluster_mp import count
D,c,M=map(int,sys.argv[1:4]); L=D-1-(D+2*c)//2
for f in sorted(glob.glob(f"rt_D{D}_c{c}_M{M}_*.npy"))[:int(sys.argv[4]) if len(sys.argv)>4 else 1]:
    A=np.load(f); Fs=A[1:].astype(int); q=np.load(f.replace("rt_","rtq_"))
    rows=[]
    for x in Fs:
        if np.linalg.matrix_rank(np.array(rows+[list(x)]))>len(rows): rows.append(list(x))
    ok,res,qq,fix,z=certify(q,rows,L); qq=np.array(qq)
    zs,mn=count(z,D,L)
    print(f,"rank",len(rows),"Krawczyk",ok,"min|Im| %.4f sep %.4f"%(np.abs(qq.imag).min(),degen(qq)),"cluster through Q:",len(zs),"rank",np.linalg.matrix_rank(np.array(zs)),flush=True)
