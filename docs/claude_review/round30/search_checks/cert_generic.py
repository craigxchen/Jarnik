import sys, pickle, numpy as np
sys.path.insert(0,'.')
from certify import certify
from basesolve import degen
D,M,c=map(int,sys.argv[1:4]); L=D-1-(D+2*c)//2
T=pickle.load(open(sys.argv[4],'rb')); R=pickle.load(open(sys.argv[5],'rb'))
done=0
for (idx,ok,best) in R:
    if not ok: continue
    E=np.array(T[idx]).T; Fs=(E[1:]!=E[0]).astype(int)
    # independent rows
    rows=[];
    for f in Fs:
        if np.linalg.matrix_rank(np.array(rows+[list(f)]))>len(rows): rows.append(list(f))
    okc,res,qq,fix,z=certify(best[2],rows,L)
    qq=np.array(qq)
    print("type",idx,"Krawczyk",okc,"res %.1e"%float(res),"min|Im| %.4f sep %.4f"%(np.abs(qq.imag).min(),degen(qq)),"fixed",[(j,str(v)) for j,v in fix],flush=True)
    if okc:
        done+=1; np.save(f"cert_D{D}_M{M}_c{c}_t{idx}.npy",qq)
        if done>=int(sys.argv[6] if len(sys.argv)>6 else 1): break
