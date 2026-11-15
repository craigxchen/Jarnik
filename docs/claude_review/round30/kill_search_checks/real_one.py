import sys, pickle, numpy as np
sys.path.insert(0, '.')
from realise import realise
t = int(sys.argv[1]); types = pickle.load(open('lane_D8_M5_T.pkl', 'rb'))
rng = np.random.default_rng(int(sys.argv[2]))
best = None
for rep in range(6):
    h, r, nf = realise(types[t], 0, 300, rng, True)
    if h and (best is None or h[0][0] > best[0]): best = h[0]
print("type", t, "best sep", best[0] if best else None)
if best: np.save('d8_type%d_q.npy' % t, best[1])
