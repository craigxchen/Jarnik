"""Search parameters (a,b,c,lam) for the k=5 norm-level family (Section 5 of ptolemy.md).

Blocks n_T(s) = s + r_T (monic linear), T subset [4], |T|>=2:
  r134=0, r124=a, r234=b, r123=c,
  r12=lam*b, r23=-lam*a, r13=((lam+1)b+(1-lam)a)/2,
  r24=-lam*c/2, r34=(2+lam)a-(1+lam)c, r14=((lam+2)b+(1-lam)c)/3.
Requirements: all roots integers = 0 mod 4 (so that s = 1 mod 4 makes every block = 1 mod 4),
all distinct, and for each small prime q some residue class of s avoids every class that
contains two or more roots (so a residue class of s gives pairwise coprime values).
"""
import itertools
from math import gcd

def roots(a, b, c, lam):
    num13 = (lam + 1)*b + (1 - lam)*a
    num24 = -lam*c
    num14 = (lam + 2)*b + (1 - lam)*c
    if num13 % 2 or num24 % 2 or num14 % 3:
        return None
    R = {'134': 0, '124': a, '234': b, '123': c, '12': lam*b, '23': -lam*a,
         '13': num13 // 2, '24': num24 // 2, '34': (2 + lam)*a - (1 + lam)*c, '14': num14 // 3}
    return R

def good(R):
    vals = list(R.values())
    if len(set(vals)) != len(vals):
        return False
    if any(v % 4 for v in vals):
        return False
    return True

def prime_factors(n):
    n = abs(n); f = set(); p = 2
    while p*p <= n:
        while n % p == 0:
            f.add(p); n //= p
        p += 1
    if n > 1:
        f.add(n)
    return f

def coprime_class(rootlist):
    """return modulus Q and residue s0 (s0 = 1 mod 4) such that all s + r are pairwise coprime for s = s0 mod Q"""
    primes = set()
    for x, y in itertools.combinations(rootlist, 2):
        primes |= prime_factors(x - y)
    primes.discard(2)
    conds = [(4, 1)]
    for q in sorted(primes):
        cnt = {}
        for r in rootlist:
            cnt[(-r) % q] = cnt.get((-r) % q, 0) + 1
        ok = [s for s in range(q) if cnt.get(s, 0) <= 1]
        if not ok:
            return None
        conds.append((q, ok[0]))
    # CRT
    Q, s0 = 1, 0
    for q, r in conds:
        # solve s = s0 mod Q, s = r mod q
        for k in range(q):
            if (s0 + Q*k) % q == r % q:
                s0 = s0 + Q*k; break
        Q *= q
    return Q, s0 % Q

def main():
  best = []
  for lam in range(-6, 7):
      if lam in (0, -1, 1, -2):
          continue
      for a in range(-40, 41, 4):
          for b in range(-40, 41, 4):
              for c in range(-40, 41, 4):
                  R = roots(a, b, c, lam)
                  if R is None or not good(R):
                      continue
                  # extra blocks: n1234 and singletons, roots chosen = 0 mod 4 distinct
                  used = set(R.values())
                  extra = []
                  x = 4
                  while len(extra) < 5:
                      for cand in (x, -x):
                          if cand not in used and len(extra) < 5:
                              extra.append(cand); used.add(cand)
                      x += 4
                  allr = list(R.values()) + extra
                  cc = coprime_class(allr)
                  if cc is None:
                      continue
                  spread = max(allr) - min(allr)
                  best.append((cc[0], spread, (a, b, c, lam), R, extra, cc))
  best.sort(key=lambda t: (t[0], t[1]))
  print(len(best), "parameter sets found")
  for row in best[:8]:
      print(row)

if __name__ == "__main__":
    main()
