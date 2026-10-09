"""Core version: is there a bijection pi : F_2^t minus 0 -> F_2^t minus 0 (rows minus one row  ->  nonconstant labels)
with (x+y).(pi(x)+pi(y)) = 1 (mod 2) for all distinct x, y?  This is exactly the condition under which the
proof of Theorem 1 applies to a Sylvester (Walsh) core.  Exhaustive backtracking."""
import sys
t = int(sys.argv[1]); M = 1 << t
els = list(range(1, M))
dot = [[bin(a & b).count("1") & 1 for b in range(M)] for a in range(M)]
pi = {}; used = set(); nodes = [0]
def bt(i):
    nodes[0] += 1
    if i == len(els):
        return True
    x = els[i]
    for v in els:
        if v in used:
            continue
        if all(dot[x ^ y][v ^ pi[y]] == 1 for y in els[:i]):
            pi[x] = v; used.add(v)
            if bt(i + 1):
                return True
            used.discard(v); del pi[x]
    return False
found = bt(0)
print("t=%d: core pairing %s  (search nodes %d)" % (t, ("FOUND " + str([pi[x] for x in els])) if found else "does NOT exist", nodes[0]))
