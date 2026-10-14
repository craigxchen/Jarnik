import sys, pickle, jdfs
ps=[5,13,17,29]; F,data=pickle.load(open("jsig_K3_37.pkl","rb"))
for sub in eval(sys.argv[1]):
    try:
        r,nodes,dt=jdfs.run(ps,3,sub,data,37,(),limit=int(sys.argv[2]))
        print(sub, "FEASIBLE" if r is not None else "INFEASIBLE (exact)", nodes, f"{dt:.1f}s", flush=True)
    except RuntimeError: print(sub,"node limit",flush=True)
