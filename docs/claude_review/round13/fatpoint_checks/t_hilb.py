import sys, time
from fplib import *
for name, L in [('syl8',2),('syl8',3),('sdc7',2),('pal7',2),('paley12',2),('pal11',2),('sdc11',2)]:
    t=time.time(); M, pats = MODELS[name](); A = box(M, pats, L)
    print(name, L, 'M', M, '#pats', len(pats), '|A|', len(A), '%.1fs'%(time.time()-t), flush=True)
    if M <= 8:
        d = M-1; H=[]; prev=0; cum=[]
        for e in range(0, 2*d+1):
            mon=[]
            for f in range(e+1): mon += reduced_monomials(d,f)
            r = rank_mod(eval_matrix(A, mon)); H.append(r-prev); prev=r
            if r == len(A): break
        print('  H =', H, flush=True)
