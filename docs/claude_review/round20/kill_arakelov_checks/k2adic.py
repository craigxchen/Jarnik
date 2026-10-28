"""2-adic local points of x^2+y^2=N (N odd, = 1 mod 4): number of classes of z=x+iy mod (1+i)^a,
versus the unit count 2^(a-1) that the K log 2 bound and the conic-only count would suggest.
Also the 2-adic local minimum sum_a E(M, classes_a) vs the binom(M,2) used by two_thirds / Theorem B."""
import math
def E(n, Q):
    if Q >= n: return 0
    b, r = divmod(n, Q); return r*math.comb(b+1, 2) + (Q-r)*math.comb(b, 2)
def classes(N, a, k):
    # z mod 2^k with |z|^2 = N mod 2^k, reduced mod (1+i)^a; (1+i)^a | w  iff  2^(a//2) | w and, if a odd, (1+i) | w/2^(a//2)
    m = 2**k; h = a//2
    seen = set()
    for x in range(m):
        for y in range(m):
            if (x*x + y*y - N) % m == 0:
                # canonical rep mod (1+i)^a: reduce mod 2^h, then if a odd keep parity class of (x' + y') of quotient
                X, Y = x % (2**h), y % (2**h)
                if a % 2 == 1:
                    qx, qy = (x - X)//(2**h), (y - Y)//(2**h)
                    seen.add((X, Y, (qx + qy) % 2))
                else:
                    seen.add((X, Y))
    return len(seen)
out = []
for N in (5, 13, 1176852625 % 2**12, 17, 29):
    row = []
    for a in range(1, 13):
        k = a//2 + 4
        row.append(classes(N, a, k))
    out.append('N mod 2^12 = %d: classes of norm-N points mod (1+i)^a, a=1..12: %s ; units 2^(a-1): %s' % (N, row, [2**(a-1) for a in range(1, 13)]))
    for M in (10, 50, 200):
        cl = row + [row[-1]*2**((a-12+1)//2) for a in range(13, 60)]
        tot = sum(E(M, c) for c in cl)
        out.append('   M=%d: 2-adic local min sum_a E(M,classes_a) = %d   vs binom(M,2) = %d   ratio %.3f' % (M, tot, math.comb(M, 2), tot/math.comb(M, 2)))
open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round20/kill_arakelov_checks/k2adic_out.txt', 'w').write('\n'.join(out)+'\n')
print('\n'.join(out))
