import sys, numpy as np
sys.path.insert(0,'/home/user/Jarnik/docs/claude_review/round8/stability_checks')
import check_fake_structure as F
q,b=11,8
H=F.normalize(F.paley1(q)); assign=[1]+list(range(1,q+1))
S=F.flipped_profile(H,b,assign)
A=(S+1)//2  # binary rows
M,r=A.shape
diffs={}
for i in range(M):
  for j in range(i+1,M):
    v=tuple(A[i]-A[j]); 
    key=max(v,tuple(-x for x in v))
    diffs.setdefault(key,[]).append((i,j))
print("M",M,"r",r,"pairs",M*(M-1)//2,"distinct +-differences",len(diffs), "max mult",max(len(x) for x in diffs.values()))
