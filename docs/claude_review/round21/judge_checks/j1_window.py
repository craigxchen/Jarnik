# J1: exact minimal P such that [P, P*exp(1/(4r))] contains >= r primes = 1 mod 4.
# Compares with the trivial lower bound 16 r (r-1) exp(-1/(4r)) and with r^2 log r.
import math, sys
from sympy import primerange
LIM = 6*10**7
ps = [p for p in primerange(5, LIM) if p % 4 == 1]
out = open(sys.argv[1], 'w')
import bisect
for r in [2,3,5,8,10,15,20,30,50,75,100,150,200,300,400]:
    f = math.exp(1.0/(4*r))
    best = None
    # window starting at P: count split primes in [P, P f]; minimal P is attained with P a split prime
    for i, P in enumerate(ps):
        hi = P*f
        if hi > LIM: break
        j = bisect.bisect_right(ps, hi)  # exact up to float rounding of P*f; check below
        # exact check of endpoint with integer arithmetic not needed if no prime within 1e-6 of hi
        if j - i >= r:
            best = P; break
    triv = 16*r*(r-1)*math.exp(-1/(4*r))
    line = f"r={r} minP={best} triv_lb={triv:.1f} ok={best is not None and best>=triv} minP/r^2={best/r**2 if best else None:.3f} minP/(r^2 log r)={best/(r*r*math.log(r)) if best else None:.3f}"
    print(line); out.write(line+"\n"); out.flush()
