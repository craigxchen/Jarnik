from fstar_q import F8
from math import log
for M in [512,1024,2048]:
    f8=F8(3,M); W=f8/M
    print(M, f8, round(W,1), round(3*M*log(M)/log(3),1), round(W/(3*M*log(M)/log(3)),3), "  (W-3Mlog3M)/M =", round((W-3*M*log(M)/log(3))/M,3), flush=True)
