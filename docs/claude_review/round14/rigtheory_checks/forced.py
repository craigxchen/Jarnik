"""Exact forced sets F_K^ex for a prime set (round 14).  Usage: python3 forced.py K qmax p1 p2 ...
Every free vector returned by the C DFS is re-verified in Python against all v in V_K."""
import sys, ctypes, os
from rlib14 import *
lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "freec.so"))
def is_free(n, nuv, K):
    r = len(nuv)
    arr = (ctypes.c_int * r)(*nuv); out = (ctypes.c_int * r)()
    res = lib.free_vector(n, r, K, arr, out)
    return list(out) if res else None
def verify(n, c, K, nuv, V):
    if n % 2 == 0 and any(c[j] % 2 != nuv[j] for j in range(len(c))): return False
    return all(sum(a * b for a, b in zip(v, c)) % n != 0 for v in V)
def forced_set(ps, K, qmax, verbose=False):
    N = 1
    for p in ps: N *= p
    V = VK_mod_sign(len(ps), K)
    F = []
    for q in primes_upto(qmax):
        if q == 2 or N % q == 0: continue
        n = (q - chi4(q)) // 4; nuv = nu(ps, q)
        c = is_free(n, nuv, K)
        if c is None: F.append(q)
        else: assert verify(n, c, K, nuv, V), (q, c)
    return F
if __name__ == "__main__":
    K = int(sys.argv[1]); qmax = int(sys.argv[2]); ps = [int(x) for x in sys.argv[3:]]
    F = forced_set(ps, K, qmax)
    print(f"ps={ps} K={K} qmax={qmax} #forced={len(F)} q*={max(F) if F else None}")
    print("forced:", F)
