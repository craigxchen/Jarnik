# local dimension (Jacobian rank) at found solutions vs the expected count of 3.2
import numpy as np, ast
from fastK import MEM, Kof
from redk import Sys
def E(v,K): return np.concatenate([np.poly(v[MEM[j]])[K].imag for j in range(4)])
for f,(mi,m0) in [('fk_50.out',(5,0)),('fk_41.out',(4,1)),('fk_42.out',(4,2))]:
    K=Kof(mi,m0)
    for line in open(f):
        if not line.startswith('SOL'): continue
        v=np.array(ast.literal_eval(line[4:].replace('np.complex128','')))
        h=1e-7; cols=[]
        for k in range(15):
            for d in (1,1j):
                e=np.zeros(15,complex); e[k]=h*d; cols.append((E(v+e,K)-E(v-e,K))/(2*h))
        J=np.array(cols).T; sv=np.linalg.svd(J,compute_uv=False); r=int(np.sum(sv>1e-7*sv[0]))
        print((mi,m0),'equations',4*len(K),'rank',r,'local dim incl. gauge',30-r,'expected',30-4*len(K))
S=Sys(4)
for line in open('k4_calib.out'):
    if 'NONDEG' in line:
        x=np.array(ast.literal_eval(line[line.index('['):].strip().replace('np.float64','')))
        h=1e-7; J=np.array([(S.residual(x+h*e)-S.residual(x-h*e))/(2*h) for e in np.eye(S.nx)]).T
        sv=np.linalg.svd(J,compute_uv=False); r=int(np.sum(sv>1e-6*sv[0]))
        print('k=4 two-branch: rank',r,'local dim (gauge fixed)',S.nx-r,'expected 4'); break
