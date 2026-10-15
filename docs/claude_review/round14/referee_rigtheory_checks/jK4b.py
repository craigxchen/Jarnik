import pickle, sys, jdfs
ps=[5,13,17,29]; K=4; data=pickle.load(open("jK4_data.pkl","rb"))
small=[11,19,23,31,37,41,43,47,53,59,61]
for q0 in [int(x) for x in sys.argv[1].split(",")]:
    sub=[3,7]+sorted(set(small+[q0]))
    try:
        r,nodes,dt=jdfs.run(ps,K,sub,data,q0,(3,7),limit=int(sys.argv[2]))
        print(f"K=4 q0={q0} G={{3,7}} subset F∩[11,61]∪{{q0}}: {'FEASIBLE(subset)' if r is not None else 'INFEASIBLE (exact)'} nodes={nodes} {dt:.0f}s",flush=True)
    except RuntimeError: print(f"K=4 q0={q0}: node limit",flush=True)
