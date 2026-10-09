import numpy as np,time
p=2147483647
rng=np.random.default_rng(0)
h=rng.random((2000,4000))*2**47
T=time.time()
for rep in range(10):
    hi=h.astype(np.int64); m=hi%p; q=hi*3
print('10 reps %.3f'%(time.time()-T))
