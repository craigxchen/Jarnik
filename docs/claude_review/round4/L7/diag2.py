import numpy as np
from newton_log import run
from newton_beta import Aps, Ts
for seed in [5,16,34]:
    ok,c,t,s,it = run(seed)
    print(seed, it, 'c', np.round(np.abs(c),4))
    vals = list(t)+list(s); names=[('t',sorted(A)) for A in Aps]+[('s',sorted(T)) for T in Ts]
    vals=np.array(vals)
    D = np.abs(vals[:,None]-vals[None,:])
    for a in range(len(vals)):
        for b in range(a+1,len(vals)):
            if D[a,b] < 1e-4: print('   close', names[a], names[b], D[a,b])
    print('  max|t|', np.max(np.abs(t)), 'max |s|', np.max(np.abs(s)))
