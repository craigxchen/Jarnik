from aux_gauss import *
import time
random.seed(1)
for primes in [[5,13,17,29,37,41,53,61,73,89], split_primes(100,200)[:11], [5,13,17,29,37,41,53,61,73,89,97,101,109]]:
    c = Circle(primes)
    t=time.time()
    for M in (3,4,5,6):
        bc = best_clusters(c, M)
        print(len(primes), "M=",M, "bestC=%.4f"%bc[0][0], [w[1] for w in bc[0][1]][:3], "time %.1f"%(time.time()-t))
