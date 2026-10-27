# prime sets used in round 14 (fixed, reproducible)
import random
from rl import split_primes
SP = split_primes(5, 400)
def sets(r):
    rng = random.Random(1400 + r)
    bal_start = SP.index(97)
    return {"first": SP[:r], "rand": sorted(rng.sample(SP, r)), "bal": SP[bal_start:bal_start + r], "last": SP[-r:]}
if __name__ == '__main__':
    for r in range(4, 9): print(r, sets(r))
