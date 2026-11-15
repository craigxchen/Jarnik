import sys, re, ast, pickle, numpy as np
sys.path.insert(0,'.'); from basesolve import res_base, degen, zero_sum_sets
from scipy.optimize import least_squares
D,L=8,3
for log,pk in (("ext8_0_9.txt","types_D8_M5_c0_sat.pkl"),("ext8_a.txt","t8_a.pkl"),("ext8_b.txt","t8_b.pkl"),("ext8_b_partial.txt","t8_b.pkl"),("ext8_a_partial.txt","t8_a.pkl")):
    T=pickle.load(open(pk,'rb'))
    for line in open(log):
        if not line.startswith("HIT"): continue
        m=re.match(r"HIT type (\d+) F \[([0-9 ]+)\] degen (\S+) (.*)",line.strip())
        idx=int(m.group(1)); F=np.array(list(map(int,m.group(2).split())),float); q=np.array(ast.literal_eval(m.group(4)))
        E=np.array(T[idx]).T; Fs=(E[1:]!=E[0]).astype(float); rt=np.linalg.matrix_rank(Fs)
        A=np.vstack([Fs,F]); z0=np.concatenate([q.real,q.imag])
        r=least_squares(res_base,z0,args=(A,L),method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=5000)
        qq=r.x[:D]+1j*r.x[D:]; S=zero_sum_sets(qq,L,1e-9)
        print(log,"type",idx,"type rank",rt,"cost %.1e"%np.linalg.norm(r.fun),"degen %.2e"%degen(qq),"cluster",len(S),"rank",np.linalg.matrix_rank(S) if len(S) else 0)
