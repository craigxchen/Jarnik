from mbar import *
import time
L8 = frozenset(range(8))
psi = {}
for i in range(8): psi[('P', i)] = Fr(1)
delta = {}
for S in all_splits(L8): delta[('D', canon(L8, S))] = Fr(1)
print(len(delta))
KD = dict(psi)
for k, v in delta.items(): KD[k] = KD.get(k, 0) - v
reps = {2: {0,1}, 3: {0,1,2}, 4: {0,1,2,3}}
t = time.time()
for s, T in reps.items():
    v = integrate(L8, [D(L8, T), KD, KD, KD, KD])
    print("D_T.(K+d)^4, |T|=", s, ":", v, time.time()-t, flush=True)
