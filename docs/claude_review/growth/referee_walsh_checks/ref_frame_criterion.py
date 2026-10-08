"""Referee check 3: Section 9.2 uses 'odd quartic support <= 4 primes' as the test for a flipped Walsh
2-frame. Definition 8.1 additionally needs (i) one flip per point (a system of distinct flip points)
and (ii) between 5 and B primes per nonzero label. Count, over all 4-point sign configurations on k
primes (columns = sign patterns of 4 points up to global sign), how often the census criterion holds but
(i) fails, and note (ii) is impossible for k < 15 varying primes."""
import itertools
pts = range(4)
chars = [(1,1,1,1),(1,-1,1,-1),(1,1,-1,-1),(1,-1,-1,1)]  # bijection F_2^2 -> points (any works for 4 pts)
def flips_ok(pattern):
    """points x such that flipping pattern at x gives a NONZERO character (up to sign)"""
    out = []
    for x in pts:
        q = list(pattern); q[x] = -q[x]
        if tuple(q) in chars[1:] or tuple(-v for v in q) in chars[1:]:
            out.append(x)
    return out
odd_patterns = [p for p in itertools.product([1,-1], repeat=4) if p[0]*p[1]*p[2]*p[3] == -1]
# census criterion holds with exactly 4 odd primes; check SDR (one flip per point and per prime)
bad = 0; tot = 0
for combo in itertools.product(odd_patterns, repeat=4):
    tot += 1
    opts = [flips_ok(p) for p in combo]
    if not any(len(set(perm)) == 4 for perm in itertools.product(*opts)):
        bad += 1
print("4 odd primes: %d of %d ordered odd-pattern 4-tuples have NO distinct flip points (criterion holds, frame fails)" % (bad, tot))
print("Def 8.1 needs >= 5 primes on each of the 3 nonzero labels of a 2-frame: >= 15 varying primes; census k ranges 8..16")
