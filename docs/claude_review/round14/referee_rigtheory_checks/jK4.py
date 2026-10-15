import pickle, jdfs, jsig
from forced import forced_set
ps=[5,13,17,29]; K=4
F,_=forced_set(ps,K,400)
data={}
for q0 in [q for q in F if 11<=q<=79]:
    o=jsig.sigs(ps,K,q0,1)
    data[q0]=dict(lev=1,allmin=jsig.minimal(list(o)),nonclass=jsig.minimal([s for s,f in o.items() if f&2]))
    r,nodes,dt=jdfs.run(ps,K,[3,7,q0],data,q0,(3,7))
    print(f"K=4 q0={q0} G={{3,7}} only q0 free: {'FEASIBLE' if r is not None else 'INFEASIBLE (exact)'} nonclass-min={len(data[q0]['nonclass'])}",flush=True)
pickle.dump(data,open("jK4_data.pkl","wb"))
