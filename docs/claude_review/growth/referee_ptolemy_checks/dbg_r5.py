import numpy as np, itertools
import rchart5 as R
rng = np.random.default_rng(3)
# replicate run but print details
import types
src = open('rchart5.py').read()
out = R.run(rng, False)
print(out[:3])
