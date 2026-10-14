from forced import forced_set
from onemod import onemod
ps=[5,13,17,29]
F,_=forced_set(ps,6,1500); print("K=6 #F(<=1500)",len(F),"max",max(F),flush=True)
for q0 in (257,263,269):
    res,w,L=onemod(ps,6,F,q0); print("K=6 q0",q0,"surv",res[0],"nonclass",res[1],"F1",res[2],"F2",res[3],"wit",w,flush=True)
