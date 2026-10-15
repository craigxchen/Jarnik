import pickle, sys, time, jdfs, jsig1
from jsig import minimal
from forced import forced_set
ps=[5,13,17,29]; K=4; data=pickle.load(open("jK4_data.pkl","rb"))
F,_=forced_set(ps,K,400)
for q in F:
    if q in (3,7) or q in data: continue
    t0=time.time(); o=jsig1.sigs1(ps,K,q)
    data[q]=dict(lev=1,allmin=minimal(list(o)),nonclass=minimal([s for s,f in o.items() if f&2]))
    print(q,"sigs",len(o),"min",len(data[q]['allmin']),f"{time.time()-t0:.0f}s",flush=True)
    pickle.dump(data,open("jK4_data.pkl","wb"))
