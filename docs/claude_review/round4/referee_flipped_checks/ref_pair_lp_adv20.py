"""Small adversarial hill-climb for Paley-20 (capacity B=2,3) of the pair-LP optimum rho (see ref_pair_lp_adversarial.py)."""
import random, sys
sys.path.insert(0, '.')
from ref_pair_lp import paley
from ref_pair_lp_adversarial import climb
rng = random.Random(20)
H = paley(19)
for B in (2, 3):
    v, a = max((climb(H, B, rng, 20) for _ in range(2)), key=lambda t: t[0])
    print(f"Paley-20 B={B}: best rho found {v} = {float(v):.4f} assignment {a}", flush=True)
