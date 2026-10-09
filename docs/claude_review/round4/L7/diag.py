import numpy as np
from newton_beta import *
for seed in [0,9,12]:
    ok,c,t,s,it = run(seed)
    print(seed, it, 'max|c|',np.max(np.abs(c)),'min|c|',np.min(np.abs(c)),'max|t|',np.max(np.abs(t)),'max|s|',np.max(np.abs(s)))
    # closest pairs among t and s
    vals = np.concatenate([t,s])
    D = np.abs(vals[:,None]-vals[None,:])+np.eye(len(vals))*1e9
    print('  min dist among t,s', D.min())
    # s close to t?
    for n,T in enumerate(Ts):
        d = np.abs(s[n]-t)
        k = np.argmin(d)
        if d[k]<1e-5: print('  s_T',sorted(T),'~ t_A',sorted(Aps[k]))
